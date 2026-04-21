#!/usr/bin/env python
'''
Module      : extract_markers
Description : Extracts markers from COMET OME-TIFF and writes comma-delimited
              output containing markers, exposure times and background channels
              compatible with input expected by mcmicro's
              background_subtraction tool.
Copyright   : (c) WEHI SODA Hub, 2025
License     : MIT
Maintainer  : Marek Cmero (@mcmero)
Portability : POSIX
'''

import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import typer
from typing import Annotated, List
from tifffile import TiffFile


def main(
    tiff_path: Annotated[Path, typer.Argument(
        help="Path to the TIFF input file."
    )],
    remove_marker: Annotated[List[str], typer.Option(
        '--remove-marker',
        '-r',
        help="List of markers to remove from the output (case sensitive)."
    )] = []
):
    """
    Extracts markers, exposure times, and background channels from COMET
    OME-TIFF
    """
    with TiffFile(tiff_path) as tiff:
        ome_metadata = tiff.ome_metadata
        imagej_metadata = getattr(tiff, 'imagej_metadata', None)

    channel_meta = []

    if ome_metadata:
        try:
            root = ET.fromstring(ome_metadata)
            ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}

            # Extract channel name and index
            channels = root.findall('.//ome:Channel', ns)
            for ch in channels:
                name = ch.attrib.get('Name')
                channel_index = ch.attrib.get('ID')
                channel_meta.append({'index': channel_index, 'marker_name': name})

            # Extract exposure times
            planes = root.findall('.//ome:Plane', ns)
            for pl in planes:
                c_index = int(pl.attrib.get('TheC'))
                exposure_time = pl.attrib.get('ExposureTime')
                if exposure_time is not None:
                    # planes are in channel order
                    channel_meta[c_index]['exposure'] = float(exposure_time)

            # Extract background channels
            channel_privs = root.findall('.//ome:ChannelPriv', ns)
            for idx, priv in enumerate(channel_privs):
                channel_id = priv.attrib.get('ID')
                for cm in channel_meta:
                    if cm['index'] == channel_id:
                        background = priv.attrib.get('FluorescenceChannel')
                        if background == cm['marker_name']:
                            channel_meta[idx]['background'] = ""
                        else:
                            channel_meta[idx]['background'] = background
                        break
        except ET.ParseError:
            channel_meta = []

    # Fallback for ImageJ metadata exports (no OME-XML)
    if not channel_meta and isinstance(imagej_metadata, dict):
        labels = imagej_metadata.get('Labels') or []
        if labels:
            channel_meta = [
                {
                    'index': str(i),
                    'marker_name': str(name),
                    'background': '',
                    'exposure': 1.0,
                }
                for i, name in enumerate(labels)
            ]
            typer.echo(
                "Warning: OME-XML metadata not found/invalid. "
                "Using ImageJ channel labels with default exposure=1.0 and blank background.",
                err=True
            )

    if not channel_meta:
        raise ValueError(
            "Could not extract channel metadata from OME-XML or ImageJ labels."
        )

    # Validate remove channel list
    markers: List[str] = [cm['marker_name'] for cm in channel_meta]
    for remove in remove_marker:
        if remove not in markers:
            raise ValueError(
                f"Channel '{remove}' not found in the OME-TIFF metadata."
            )

    # Ensure required fields exist in fallback mode
    for cm in channel_meta:
        if 'background' not in cm or cm['background'] is None:
            cm['background'] = ''
        if 'exposure' not in cm or cm['exposure'] is None:
            cm['exposure'] = 1.0

    # Format output and write to stdout
    df = pd.DataFrame(channel_meta)
    df = df[['marker_name', 'background', 'exposure']]
    df['exposure'] = df['exposure'].astype(float)

    # Handle removal of channels
    markers_to_remove: List[str] = \
        ['TRUE' if marker in remove_marker else ''
            for marker in df['marker_name']]
    df['remove'] = markers_to_remove

    df.to_csv(sys.stdout, index=False)


if __name__ == "__main__":
    typer.run(main)
