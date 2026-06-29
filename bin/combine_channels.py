#!/usr/bin/env python
'''
Module      : combine_channels
Description : Takes an OME-TIFF (including OPAL OME/QPTIFF) containing N
              channels and returns a 2 channel tiff
              containing a nuclear channel and membrane channel that is componsed
              of one or more channels from the input tiff, usign either the product
              or max of the intensities.
Copyright   : (c) WEHI SODA Hub, 2025
License     : MIT
Maintainer  : Marek Cmero (@mcmero)
Portability : POSIX
'''
import json
import re
import sys
import xml.etree.ElementTree as ET
from enum import Enum
from pathlib import Path
from typing import Annotated, Any, Dict, List, Optional, Tuple

import typer
import numpy as np
from tifffile import TiffFile, imwrite
from xarray import DataArray


class CombineMethod(str, Enum):
    PROD = "prod"
    MAX = "max"


def _local_tag_name(tag: str) -> str:
    """
    Return XML tag name without namespace.
    """
    return tag.split('}', 1)[-1]


def get_pixels_tag(xml_str: str) -> ET.Element:
    """
    Parses the OME-XML string and returns the Pixels tag.
    """
    root = ET.fromstring(xml_str)

    image_tag: Optional[ET.Element] = None
    for elem in root.iter():
        if _local_tag_name(elem.tag) == 'Image':
            image_tag = elem
            break

    if image_tag is None:
        raise ValueError("No Image tag found in the XML.")

    for child in image_tag:
        if _local_tag_name(child.tag) == 'Pixels':
            return child

    raise ValueError("No Pixels tag found in the Image tag.")


def ome_extract_channel_names(xml_str: str) -> List[str]:
    """
    Extracts all channel 'Name' attributes from OME-XML in the order they
    appear. Returns a list of channel names by index.
    """
    pixels_tag = get_pixels_tag(xml_str)

    channel_tags = [
        ch for ch in list(pixels_tag)
        if _local_tag_name(ch.tag) == 'Channel'
    ]
    channel_names = [ch.get('Name', ch.get('ID', 'Unknown')) for ch in channel_tags]
    return channel_names


def ome_extract_pixels_metadata(xml_str: str) -> Dict[str, Any]:
    """
    Extract commonly used metadata fields from OME-XML Pixels tag.
    """
    pixels_tag = get_pixels_tag(xml_str)

    fields = [
        "DimensionOrder", "Type", "SizeX", "SizeY", "SizeZ", "SizeC", "SizeT",
        "PhysicalSizeX", "PhysicalSizeY", "PhysicalSizeXUnit", "PhysicalSizeYUnit",
    ]
    return {f: pixels_tag.get(f) for f in fields if pixels_tag.get(f) is not None}


def json_extract_channel_names(pages) -> List[str]:
    """
    Extracts channel names from MIBI TIFF pages where each page has a JSON
    description containing channel metadata.
    """
    channel_names: List[str] = []
    for page in pages:
        desc = json.loads(page.description)
        channel_names.append(desc["channel.target"])

    return channel_names


def imagej_extract_channel_names(tiff: TiffFile) -> List[str]:
    """
    Extract channel names from ImageJ metadata when available.
    """
    imagej_metadata = getattr(tiff, "imagej_metadata", None)
    if not imagej_metadata or not isinstance(imagej_metadata, dict):
        return []

    labels = imagej_metadata.get("Labels")
    if not labels:
        return []

    return [str(label) for label in labels]


def detect_channel_names_and_metadata(tiff: TiffFile) -> Tuple[List[str], Dict[str, Any], str]:
    """
    Detect channel names and metadata from supported TIFF metadata formats.

    Returns:
        (channel_names, metadata_attrs, metadata_source)
    """
    channel_names: List[str] = []
    metadata_attrs: Dict[str, Any] = {}
    metadata_source = "unknown"

    # MIBI JSON metadata in per-page ImageDescription
    first_page = tiff.pages[0]
    try:
        json.loads(first_page.description)
        channel_names = json_extract_channel_names(tiff.pages)
        metadata_source = "mibi-json"
    except (json.JSONDecodeError, TypeError):
        pass

    # OME metadata (covers standard OME-TIFF and OPAL OME/QPTIFF exports)
    if not channel_names:
        ome_xml = tiff.ome_metadata
        if not ome_xml and isinstance(first_page.description, str):
            description = first_page.description.strip()
            if description.startswith("<"):
                ome_xml = first_page.description

        if ome_xml:
            try:
                channel_names = ome_extract_channel_names(ome_xml)
                metadata_attrs.update(ome_extract_pixels_metadata(ome_xml))
                metadata_source = "ome-xml"
            except (ET.ParseError, ValueError, TypeError):
                pass

    # ImageJ metadata (e.g. OPAL exports without embedded OME-XML)
    if not channel_names:
        channel_names = imagej_extract_channel_names(tiff)
        if channel_names:
            metadata_source = "imagej"

    return channel_names, metadata_attrs, metadata_source


def tiff_to_xarray(tiffPath: Path) -> DataArray:
    """
    Takes a TIFF and converts it to an xarray with relevant axis,
    coordinate and metadata attached. Supports MIBI TIFF, OME-TIFF,
    and OPAL OME/QPTIFF exports.
    Uses memory mapping to avoid loading the entire image into memory.
    """
    channel_names: List[str] = []
    attrs: Dict[str, Any] = {}

    with TiffFile(tiffPath) as tiff:
        channel_names, extracted_attrs, metadata_source = detect_channel_names_and_metadata(tiff)
        attrs.update(extracted_attrs)
        attrs["metadata_source"] = metadata_source

        if len(tiff.pages) > 1:
            # Stack pages using memory mapping
            arrays = [page.asarray() for page in tiff.pages]
            data = np.stack(arrays, axis=0)
        else:
            data = tiff.asarray()
            # Single-page TIFFs with interleaved channels have shape (Y, X, C).
            # Transpose to (C, Y, X) if the last axis matches known channel names,
            # or if it looks like a channel axis (much smaller than spatial dims)
            # when metadata is absent.
            if data.ndim == 3:
                last_matches_names = bool(channel_names) and data.shape[2] == len(channel_names)
                looks_like_channel_last = (
                    not channel_names and data.shape[2] < data.shape[0] and data.shape[2] < data.shape[1]
                )
                if last_matches_names or looks_like_channel_last:
                    data = np.transpose(data, (2, 0, 1))

        if not channel_names and data.ndim == 3:
            channel_names = [f"Channel_{i}" for i in range(data.shape[0])]

        if not channel_names:
            raise ValueError(
                "Could not detect channel names from TIFF metadata. "
                "Please provide an image with OME or MIBI metadata."
            )

        return DataArray(data=data, dims=["C", "Y", "X"],
                         coords={"C": channel_names}, attrs=attrs)


def combine_channels(array: DataArray, channels: List[str], combined_name: str,
                     combine_method: CombineMethod) -> DataArray:
    """
    Combines multiple channels into a single channel using the specified method
    (prod or max). Adds the combined channel to the array.
    """

    if len(channels) == 1:
        return array

    # Select the specified channels
    selected_data = array.sel(C=channels).values

    if combine_method == CombineMethod.MAX:
        combined_data = np.max(selected_data, axis=0, keepdims=True)
    elif combine_method == CombineMethod.PROD:
        # Convert to uint64 to avoid possible overflow
        selected_data = selected_data.astype(np.uint64)
        combined_data = np.prod(selected_data, axis=0, keepdims=True)

    # Convert back to uint16, scaling if necessary
    max_val: int = np.iinfo(np.uint16).max
    if np.max(combined_data) > max_val:
        scale_factor: float = (np.iinfo(np.uint16).max - 1) / max_val
        combined_data = np.clip(combined_data * scale_factor, 0,
                                np.iinfo(np.uint16).max).astype(np.uint16)
    else:
        combined_data = combined_data.astype(np.uint16)

    # Create new array with combined channel
    new_data = np.concatenate([array.values, combined_data], axis=0)
    new_coords = list(array.coords["C"].values) + [combined_name]

    # Delete intermediate arrays to free memory
    del selected_data, combined_data

    return DataArray(data=new_data, dims=["C", "Y", "X"],
                     coords={"C": new_coords}, attrs=array.attrs)


def create_ome_xml(width: int, height: int,
                   num_channels: int, channel_names: List[str],
                   pixel_size_microns: Optional[float] = None) -> str:
    """
    Create minimal OME-XML metadata from scratch with the given channel names.
    Used when the source image has no OME-XML (e.g. MIBI TIFF with JSON metadata).
    """
    ome_ns = "http://www.openmicroscopy.org/Schemas/OME/2016-06"
    xsi_ns = "http://www.w3.org/2001/XMLSchema-instance"
    schema_loc = (
        "http://www.openmicroscopy.org/Schemas/OME/2016-06 "
        "http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd"
    )
    ET.register_namespace("", ome_ns)
    ET.register_namespace("xsi", xsi_ns)

    # Pre-build Clark-notation tag names to avoid {{{ patterns that trip nf-core linting
    tag_ome = "{" + ome_ns + "}"
    tag_xsi = "{" + xsi_ns + "}"

    root = ET.Element(tag_ome + "OME")
    root.set(tag_xsi + "schemaLocation", schema_loc)

    image = ET.SubElement(root, tag_ome + "Image", {"ID": "Image:0", "Name": "combined"})
    pixels_attrs = {
        "ID": "Pixels:0",
        "DimensionOrder": "XYZCT",
        "Type": "uint16",
        "SizeX": str(width),
        "SizeY": str(height),
        "SizeZ": "1",
        "SizeC": str(num_channels),
        "SizeT": "1",
    }
    if pixel_size_microns is not None:
        pixels_attrs["PhysicalSizeX"] = str(pixel_size_microns)
        pixels_attrs["PhysicalSizeY"] = str(pixel_size_microns)
        pixels_attrs["PhysicalSizeXUnit"] = "µm"
        pixels_attrs["PhysicalSizeYUnit"] = "µm"
    pixels = ET.SubElement(image, tag_ome + "Pixels", pixels_attrs)
    for c, name in enumerate(channel_names):
        ET.SubElement(pixels, tag_ome + "Channel", {
            "ID": f"Channel:0:{c}",
            "Name": name,
            "SamplesPerPixel": "1",
        })
    tiff_data = ET.SubElement(pixels, tag_ome + "TiffData")
    tiff_data.set("PlaneCount", str(num_channels))

    xml_str = '<?xml version="1.0" encoding="UTF-8"?>\n' + \
        ET.tostring(root, encoding='unicode')
    return xml_str


def update_ome_xml(original_xml: str, width: int, height: int,
                   num_channels: int, channel_names: List[str]) -> str:
    """
    Update existing OME-XML metadata with new channel information.
    """
    try:
        if isinstance(original_xml, bytes):
            original_xml = original_xml.decode('utf-8', errors='replace')

        original_xml = re.sub(r'<\?xml[^>]+\?>', '<?xml version="1.0"?>', original_xml)
        root = ET.fromstring(original_xml)

        namespace = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}
        pixels = root.find('.//ome:Pixels', namespace)

        if pixels is not None:
            # Update dimensions
            pixels.set('SizeX', str(width))
            pixels.set('SizeY', str(height))
            pixels.set('SizeC', str(num_channels))

            # Update/set PlaneCount to match channel number
            tiff_data_elements = pixels.findall('ome:TiffData', namespace)
            if tiff_data_elements:
                for tiff_data in tiff_data_elements:
                    tiff_data.set('PlaneCount', str(num_channels))
            else:
                tiff_data = ET.SubElement(pixels, 'TiffData')
                tiff_data.set('PlaneCount', str(num_channels))

            for channel in pixels.findall('ome:Channel', namespace):
                pixels.remove(channel)

            for c in range(num_channels):
                channel_id = f"Channel:{c}"
                channel_attrs = {
                    'ID': channel_id,
                    'Name': channel_names[c],
                    'SamplesPerPixel': '1'
                }
                ET.SubElement(pixels, 'Channel', channel_attrs)

        xml_str = '<?xml version="1.0" encoding="UTF-8"?>\n' + \
            ET.tostring(root, encoding='utf-8').decode('utf-8')

        # Fix ns0 namespace prefixes if they exists
        xml_str = re.sub(r'<ns0:', '<', xml_str)
        xml_str = re.sub(r'</ns0:', '</', xml_str)
        xml_str = re.sub(r'xmlns:ns0', 'xmlns', xml_str)

        return xml_str

    except ET.ParseError as e:
        raise ValueError(f"Failed to parse OME-XML metadata: {e}")


def main(
    tiff: Annotated[Path, typer.Argument(
        help="Path to the TIFF input file."
    )],
    nuclear_channel: Annotated[str, typer.Option(
        help="Name of the nuclear channel."
    )],
    membrane_channel: Annotated[List[str], typer.Option(
        help="Name(s) of the membrane channels (can be repeated)"
             "Ensure that channels with spaces are quoted.")
    ],
    combine_method: Annotated[CombineMethod, typer.Option(
        help="Method to use for combining channels (prod or max).")
    ] = CombineMethod.MAX,
    pixel_size_microns: Annotated[Optional[float], typer.Option(
        help="Pixel size in microns. Written to OME-XML PhysicalSizeX/Y when "
             "the source image has no calibration metadata.")
    ] = None,
):
    full_array = tiff_to_xarray(tiff)

    # Combine membrane channels if needed
    if len(membrane_channel) > 1:
        combined_membrane_channel = "combined_membrane"
        full_array = combine_channels(full_array, membrane_channel,
                                      combined_membrane_channel,
                                      CombineMethod(combine_method))
        final_channels = [nuclear_channel, combined_membrane_channel]
    else:
        final_channels = [nuclear_channel, membrane_channel[0]]

    # Report what was combined so it stays visible in the process log
    # (downstream segmenters only see the single combined membrane channel).
    typer.echo(
        f"[combine_channels] nuclear='{nuclear_channel}'; "
        f"membrane={membrane_channel} -> '{final_channels[1]}' "
        f"(method={CombineMethod(combine_method).value}, "
        f"n={len(membrane_channel)})",
        err=True,
    )

    # Extract final channels and convert to output format
    output_array = full_array.sel(C=final_channels).values.astype(np.uint16)

    # Free the full_array from memory
    del full_array

    with TiffFile(tiff) as tif:
        ome_metadata = tif.ome_metadata or tif.pages[0].description

    # Update OME-XML metadata
    c, height, width = output_array.shape
    try:
        updated_metadata = update_ome_xml(ome_metadata, width, height, c,
                                          final_channels)
        updated_metadata = updated_metadata.encode('utf-8')
    except ValueError as e:
        typer.echo(f"Warning: {e} Creating OME-XML from scratch.",
                   err=True)
        updated_metadata = create_ome_xml(width, height, c,
                                          final_channels,
                                          pixel_size_microns).encode('utf-8')

    imwrite(sys.stdout.buffer, output_array,
            photometric='minisblack',
            metadata={'axes': 'CYX'},
            description=updated_metadata)


if __name__ == "__main__":
    typer.run(main)
