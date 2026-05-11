#!/usr/bin/env python3
"""
Python rewrite of the cellmeasurement Groovy app from:
https://github.com/WEHI-SODA-Hub/cellmeasurement

This script measures cell and nucleus compartments from labeled masks and a
multi-channel TIFF image, then exports a GeoJSON FeatureCollection.
"""

from __future__ import annotations

import argparse
import gc
import gzip
import json
import math
import os
import sys
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import tifffile
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from shapely.geometry import Polygon, mapping, shape
from shapely.ops import unary_union
from skimage.measure import block_reduce, find_contours, regionprops
from skimage.morphology import binary_erosion, disk
from skimage.segmentation import watershed

# Pre-computed 3x3 disk structuring element used throughout for single-pixel
# morphological erosion/dilation operations.  Frozen to prevent accidental mutation.
_DISK_1 = disk(1)
_DISK_1.flags.writeable = False

# ---------------------------------------------------------------------------
# Module-level globals for sharing large arrays with worker processes via
# fork() copy-on-write.  Set in main() before ProcessPoolExecutor is created
# so that forked workers inherit the arrays without pickling per-task copies.
# Workers only READ from these arrays, so pages are never actually copied.
# ---------------------------------------------------------------------------
_GLOBAL_IMG: Optional[np.ndarray] = None    # full multi-channel image (original dtype)
_GLOBAL_CELL: Optional[np.ndarray] = None   # whole-cell label mask
_GLOBAL_NUC: Optional[np.ndarray] = None    # nuclear label mask (may be None)
_GLOBAL_SKIP_NUC: bool = False


@dataclass
class CellRecord:
    """Association record linking a unified cell ID to its source mask labels.

    After ``match_cells()`` pairs nuclear and whole-cell segmentation labels,
    each matched (or synthesized) cell gets a unique ``cell_id`` used in the
    output label arrays.  ``cell_label`` and ``nucleus_label`` refer back to
    the *original* label values in the input whole-cell and nuclear masks
    respectively, allowing provenance tracking.  Either may be ``None`` when
    the cell was synthesized from only one mask.
    """
    cell_id: int
    cell_label: Optional[int]
    nucleus_label: Optional[int]


def parse_args() -> argparse.Namespace:
    """Define and parse all command-line arguments.

    Returns a namespace with validated arguments controlling:

    * **Input/output paths** — nuclear mask, whole-cell mask, multi-channel
      TIFF image, and GeoJSON output file.
    * **Spatial parameters** — pixel size, downsample factor, distance
      thresholds for nucleus-to-cell matching and synthetic boundary
      estimation.
    * **Measurement options** — percentile computation, erosion/expansion
      ring steps, environment dilation distance, and neighbourhood
      aggregation (k-nearest neighbours).
    * **Geometry options** — ROI simplification toggle and Douglas-Peucker
      tolerance.
    * **Performance** — thread count, tile size/overlap (compatibility only).
    * **Output variants** — pretty-printed JSON, optional rasterized label
      mask TIFF, and ``--skip-nuclear-mask`` mode.
    """
    p = argparse.ArgumentParser(description="Extract cell measurements from masks and image TIFF.")
    p.add_argument("-n", "--nuclear-mask", required=True, help="Nuclear segmentation mask TIFF")
    p.add_argument("-w", "--whole-cell-mask", required=True, help="Whole-cell segmentation mask TIFF")
    p.add_argument("-f", "--tiff-file", required=True, help="Multi-channel TIFF image")
    p.add_argument("-o", "--output-file", required=True, help="Output GeoJSON path")
    p.add_argument("-d", "--downsample-factor", type=float, default=1.0)
    p.add_argument("-p", "--pixel-size-microns", type=float, default=0.5)
    p.add_argument("--skip-measurements", action="store_true")
    p.add_argument("--simplify-rois", action="store_true",
                   help="Simplify ROI geometry with Douglas-Peucker. Enabled by default; use --no-simplify-rois to disable.")
    p.add_argument("--no-simplify-rois", dest="simplify_rois", action="store_false")
    p.set_defaults(simplify_rois=True)
    p.add_argument("--tolerance", type=float, default=0.5,
                   help="Simplification tolerance in pixels (default 0.5). Lower values preserve more shape detail.")
    p.add_argument("--percentiles", default="")
    p.add_argument("--erosion-steps", action="store_true",
                   help="Measure intensity in 5 equal-area erosion bins working inward from the cell boundary. "
                        "Disabled by default. "
                        "NOTE: unlike the original Groovy cellmeasurement app, this flag does not accept pixel-depth "
                        "values (e.g. 4,7,11,14,18). The Python implementation always uses 5 equal-area bins; "
                        "the output measurement keys are named ErosionBin_1 through ErosionBin_5, not Eroded_Npx.")
    p.add_argument("--expansion-steps", action="store_true",
                   help="Measure intensity in 5 equal-area annular bins dilated 20 µm outward from cell boundary "
                        "(distance computed from --pixel-size-microns). Disabled by default.")
    p.add_argument("-i", "--dist-threshold", type=float, default=10.0)
    p.add_argument("-e", "--estimate-cell-boundary-dist", type=float, default=3.0)
    p.add_argument("-t", "--threads", type=int, default=1)
    p.add_argument("--tile-size", type=int, default=2048)
    p.add_argument("--tile-overlap", type=int, default=200)
    p.add_argument("--pretty-json", action="store_true", help="Write indented GeoJSON output")
    p.add_argument("--gzip", action="store_true",
                   help="Gzip-compress the output GeoJSON file (appends .gz to filename if needed)")
    p.add_argument("--output-mask", default="",
                   help="Write a rasterized label mask TIFF from the final cell geometries")
    p.add_argument("--skip-nuclear-mask", action="store_true",
                   help="Ignore the nuclear mask entirely. ROIs are generated from the whole-cell mask only. "
                        "Compartmental measurements (Nucleus, Cytoplasm, Membrane) are skipped. "
                        "Cell, erosion, and expansion measurements are still produced.")
    p.add_argument("--neighbors", type=int, default=0,
                   help="Number of nearest neighbors for neighborhood feature aggregation (0 = disabled). "
                        "Computes max and mean of every numeric measurement across each cell's k closest neighbours.")
    p.add_argument("--environment-expansion", action="store_true",
                   help="Measure a pericellular 'Environment' compartment by dilating the cell mask "
                        "outward by 20 µm (converted to pixels via --pixel-size-microns). "
                        "Disabled by default.")
    return p.parse_args()


def tile_flags_explicit(argv: Sequence[str]) -> bool:
    """Check whether the user explicitly passed ``--tile-size`` or ``--tile-overlap``.

    These flags are accepted for CLI compatibility with the original Groovy
    implementation but are **not used** in this Python version.  When detected
    a warning is printed so the user knows the flags are no-ops.

    Parameters
    ----------
    argv : Sequence[str]
        The raw command-line arguments (typically ``sys.argv[1:]``).

    Returns
    -------
    bool
        True if either flag appears in *argv*.
    """
    return any(
        a == "--tile-size"
        or a.startswith("--tile-size=")
        or a == "--tile-overlap"
        or a.startswith("--tile-overlap=")
        for a in argv
    )


def load_label_mask(path: str) -> np.ndarray:
    """Read a 2-D integer label mask from a TIFF file.

    The mask is expected to be a single-plane image where each pixel's value
    is the integer label of the segmented object (0 = background).  The array
    is cast to ``int64`` so that downstream arithmetic (differences, products)
    never overflows, regardless of the original bit depth (commonly uint16 or
    uint32 from segmentation tools).

    Parameters
    ----------
    path : str
        Filesystem path to the TIFF label mask.

    Returns
    -------
    np.ndarray
        2-D ``int64`` array of shape (H, W).

    Raises
    ------
    ValueError
        If the image is not exactly 2-D.
    """
    arr = tifffile.imread(path)
    if arr.ndim != 2:
        raise ValueError(f"Mask must be 2D label image, got shape={arr.shape} for {path}")
    # Normalize to signed integer labels to keep downstream ops consistent.
    return arr.astype(np.int64, copy=False)


def load_image(path: str) -> Tuple[np.ndarray, List[str]]:
    """Load a multi-channel TIFF image and extract per-channel names.

    Supports several multiplexed imaging metadata conventions:

    1. **OME-XML** ``Channel/@Name`` — standard for OME-TIFF; also works for
       OPAL QPTIFF and COMET images that embed OME-XML in the first page's
       ``ImageDescription`` tag.
    2. **MIBI JSON** — per-page ``ImageDescription`` containing JSON with a
       ``channel.target`` field (IONpath MIBIscope convention).
    3. **ImageJ Labels** — ``imagej_metadata['Labels']`` list.

    If none of these strategies yield the correct number of names, fallback
    names ``"Channel 1"``, ``"Channel 2"``, … are generated.

    The returned array is always shaped ``(C, H, W)`` and cast to
    ``float32`` for downstream intensity arithmetic.

    Parameters
    ----------
    path : str
        Filesystem path to the multi-channel TIFF.

    Returns
    -------
    image : np.ndarray
        ``float32`` array of shape ``(C, H, W)``.
    ch_names : list of str
        Human-readable channel names, length ``C``.
    """
    # Use the same loading strategy as cellsam_segment.py to ensure channel
    # detection stays consistent across modules.  Iterating over ``tif.pages``
    # explicitly avoids ``tifffile.imread`` collapsing a multi-page TIFF
    # (e.g. an OPAL QPTIFF or BACKSUB output) into a single 2-D plane.
    import json as _json

    ch_names: List[str] = []

    with tifffile.TiffFile(path) as tf:
        first_page = tf.pages[0] if tf.pages else None
        is_mibi = False

        # Detect MIBI-style per-page JSON metadata.
        if first_page is not None:
            try:
                first_desc = _json.loads(first_page.description)
                is_mibi = "channel.target" in first_desc
            except (ValueError, TypeError):
                pass

        if is_mibi:
            pages = []
            for page in tf.pages:
                desc = _json.loads(page.description)
                ch_names.append(str(desc.get("channel.target", "")))
                pages.append(page.asarray())
            img = np.stack(pages, axis=0)
        elif len(tf.pages) > 1:
            # Multi-page TIFF where each page is one channel (OME-TIFF, OPAL
            # QPTIFF, COMET, BACKSUB output, etc.).
            img = np.stack([page.asarray() for page in tf.pages], axis=0)
        else:
            # Single-page TIFF — could be (H, W) or interleaved (H, W, C).
            img = tf.asarray()

        # Strategy 1: OME-XML Channel/@Name (OME-TIFF, OPAL QPTIFF, COMET)
        if not ch_names:
            try:
                import xml.etree.ElementTree as ET

                ome = tf.ome_metadata
                # OPAL QPTIFF stores OME-XML in first page ImageDescription
                if not ome and first_page is not None:
                    first_desc = first_page.description
                    if isinstance(first_desc, str) and first_desc.strip().startswith("<"):
                        ome = first_desc
                if ome:
                    root = ET.fromstring(ome)
                    ns = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
                    channels = root.findall(".//ome:Channel", ns)
                    if channels:
                        ch_names = [
                            ch.get("Name") or ch.get("ID") or ""
                            for ch in channels
                        ]
            except Exception as e:
                print(f"Warning: failed to parse OME metadata for channel names ({e})")
                ch_names = []

        # Strategy 2: ImageJ metadata Labels
        if not ch_names:
            try:
                ij = tf.imagej_metadata
                if ij and "Labels" in ij:
                    ch_names = [str(lbl) for lbl in ij["Labels"]]
            except Exception:
                pass

        n_pages = len(tf.pages) if tf.pages else 0

    # --- Normalize image shape to (C, H, W) regardless of input layout ---
    if img.ndim == 2:
        # Single-channel greyscale: promote to (1, H, W)
        img = img[np.newaxis, ...]
    elif img.ndim == 3:
        # If we stacked pages above, the array is already (C, H, W).
        # Only single-page TIFFs can land here as interleaved (H, W, C).
        if n_pages > 1 and img.shape[0] == n_pages:
            pass  # already (C, H, W) — each page is one channel
        elif n_pages == 1 or (n_pages > 1 and img.shape[2] == n_pages):
            img = np.transpose(img, (2, 0, 1))  # (H, W, C) -> (C, H, W)
        elif img.shape[0] < img.shape[1] and img.shape[0] < img.shape[2]:
            pass  # heuristic: first dim is smaller than spatial dims -> (C, H, W)
        elif img.shape[2] < img.shape[0] and img.shape[2] < img.shape[1]:
            img = np.transpose(img, (2, 0, 1))  # heuristic: last dim smallest -> (H, W, C)
        else:
            raise ValueError(
                f"Unsupported 3D image layout: shape={img.shape}, n_pages={n_pages}. "
                "Cannot determine whether channels are first or last."
            )
    else:
        raise ValueError(f"Unsupported image dimensions: {img.shape}")

    if not ch_names or len(ch_names) != img.shape[0]:
        print(f"Warning: found {len(ch_names)} channel names but image has {img.shape[0]} channels; using fallback names")
        ch_names = [f"Channel {i + 1}" for i in range(img.shape[0])]
    else:
        print(f"Detected channel names: {ch_names}")

    # Return in original dtype (typically uint16).  Casting to float32 here
    # would double the memory footprint (~70 GB for a 34-channel whole-slide
    # image).  Conversion to float32 is deferred to per-crop workers.
    return img, ch_names


def maybe_downsample(image_cyx: np.ndarray, nuc: np.ndarray, whole: np.ndarray, ds: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Optionally downsample the image and label masks by an integer factor.

    When ``ds > 1``, the image intensities are reduced using **block mean**
    (to preserve average signal) while label masks use **block max** (to
    avoid creating spurious zero-gaps between adjacent labels).  The arrays
    are first cropped to dimensions divisible by the step so that
    ``block_reduce`` does not zero-pad edge blocks and introduce bias.

    Parameters
    ----------
    image_cyx : np.ndarray
        Multi-channel image, shape ``(C, H, W)``.
    nuc : np.ndarray
        Nuclear label mask, shape ``(H, W)``.
    whole : np.ndarray
        Whole-cell label mask, shape ``(H, W)``.
    ds : float
        Downsample factor.  Rounded to the nearest integer; values < 2
        result in no downsampling.

    Returns
    -------
    tuple of (image_ds, nuc_ds, whole_ds)
        Downsampled arrays.  ``image_ds`` is ``float32``; masks are ``int64``.
    """
    if ds <= 1.0:
        return image_cyx, nuc, whole
    step = int(round(ds))
    if step < 2:
        return image_cyx, nuc, whole
    # Crop to divisible shape first to avoid block_reduce zero-padding bias at edges.
    c, h, w = image_cyx.shape
    h2 = (h // step) * step
    w2 = (w // step) * step
    if h2 == 0 or w2 == 0:
        return image_cyx, nuc, whole

    image_crop = image_cyx[:, :h2, :w2]
    nuc_crop = nuc[:h2, :w2]
    whole_crop = whole[:h2, :w2]

    # Use reshape+mean for image intensities (faster for this regular block case)
    # and block max for label masks.
    image_ds = image_crop.reshape(c, h2 // step, step, w2 // step, step).mean(axis=(2, 4))
    nuc_ds = block_reduce(nuc_crop, block_size=(step, step), func=np.max)
    whole_ds = block_reduce(whole_crop, block_size=(step, step), func=np.max)
    return image_ds.astype(np.float32, copy=False), nuc_ds.astype(np.int64, copy=False), whole_ds.astype(np.int64, copy=False)


def label_props_dict(label_img: np.ndarray) -> Dict[int, Dict[str, Any]]:
    """Build a dictionary of region properties keyed by label value.

    Uses ``skimage.measure.regionprops`` to extract the centroid (row, col)
    and bounding box (minr, minc, maxr, maxc) of each labelled region.
    These are the only properties needed by the matching and task-generation
    stages, so we extract them once upfront rather than calling
    ``regionprops`` repeatedly.

    Parameters
    ----------
    label_img : np.ndarray
        2-D integer label image (0 = background).

    Returns
    -------
    dict
        ``{label: {"centroid": (row, col), "bbox": (minr, minc, maxr, maxc)}}``
    """
    out: Dict[int, Dict[str, Any]] = {}
    for r in regionprops(label_img):
        out[int(r.label)] = {
            "centroid": (float(r.centroid[0]), float(r.centroid[1])),
            "bbox": (int(r.bbox[0]), int(r.bbox[1]), int(r.bbox[2]), int(r.bbox[3])),
        }
    return out


def match_cells(
    nuc: np.ndarray, whole: np.ndarray, dist_threshold: float, estimate_dist: float
) -> Tuple[np.ndarray, np.ndarray, List[CellRecord], Dict[str, int], Dict[int, Tuple[slice, slice]]]:
    """Match nuclear labels to whole-cell labels and produce unified cell arrays.

    This is the core cell-association step that links nuclei to their
    enclosing cytoplasm.  It runs in two passes:

    **Pass 1 — centroid matching:**
    For each nucleus, find the nearest whole-cell centroid within
    ``dist_threshold`` pixels using a ``cKDTree``.  If a unique match is
    found, the nucleus and whole-cell pixels are painted into the output
    arrays with a shared ``cell_id``.  Pixels already claimed by an
    earlier cell are excluded (first-come-first-served) to avoid label
    collisions in the output.

    **Pass 2 — watershed synthesis:**
    Nuclei that did not match any whole-cell label get synthetic cell
    boundaries via watershed expansion.  Each unmatched nucleus is used
    as a seed; the watershed floods outward up to ``estimate_dist``
    pixels into unclaimed territory, producing Voronoi-like partitions
    that prevent adjacent unmatched nuclei from overlapping.

    Parameters
    ----------
    nuc : np.ndarray
        Nuclear label mask (H, W), int64.
    whole : np.ndarray
        Whole-cell label mask (H, W), int64.
    dist_threshold : float
        Maximum centroid distance (pixels) for a nucleus–cell match.
    estimate_dist : float
        Dilation radius (pixels) for synthesizing cell boundaries around
        unmatched nuclei.

    Returns
    -------
    out_cell : np.ndarray
        Unified cell label mask (H, W) with sequential ``cell_id`` values.
    out_nuc : np.ndarray
        Corresponding nuclear label mask (H, W), using the same
        ``cell_id`` values as ``out_cell``.
    records : list of CellRecord
        One record per cell, linking ``cell_id`` to the original mask labels.
    stats : dict
        Matching summary counts (nuclei, whole cells, matched, unmatched,
        dropped).
    bbox_map : dict
        ``{cell_id: (row_slice, col_slice)}`` tight bounding box for each
        cell, padded by 1 pixel on each side for safe contour extraction.
    """
    whole_props = label_props_dict(whole)
    nuc_props = label_props_dict(nuc)

    records: List[CellRecord] = []
    bbox_map: Dict[int, Tuple[slice, slice]] = {}
    out_cell = np.zeros_like(whole, dtype=np.int64)
    out_nuc = np.zeros_like(nuc, dtype=np.int64)

    # Build cKDTree from whole-cell centroids for fast nearest-neighbour lookup
    # when matching each nucleus to its enclosing cell body.
    whole_labels = sorted(whole_props.keys())
    whole_pts = np.array([whole_props[l]["centroid"] for l in whole_labels], dtype=np.float64)
    tree = cKDTree(whole_pts) if len(whole_pts) else None

    next_id = 1  # sequential cell ID counter for the output arrays
    used_whole: set[int] = set()  # whole-cell labels already claimed by a nucleus
    dropped_synth_cells = 0  # count of cells dropped due to complete occlusion

    # Track which nuclei need synthesized boundaries for the watershed pass.
    unmatched_nuclei: List[Tuple[int, int]] = []  # (nuc_label, assigned_cell_id)

    # --- Pass 1: match nuclei to whole-cell masks; defer unmatched nuclei ---
    for nlab in sorted(nuc_props.keys()):
        nr, nc = nuc_props[nlab]["centroid"]
        nminr, nminc, nmaxr, nmaxc = nuc_props[nlab]["bbox"]
        matched_whole = None

        if tree is not None:
            dist, idx = tree.query([nr, nc], k=1, distance_upper_bound=dist_threshold)
            if np.isfinite(dist) and idx < len(whole_labels):
                candidate = whole_labels[int(idx)]
                if candidate not in used_whole:
                    matched_whole = candidate

        if matched_whole is not None:
            cminr, cminc, cmaxr, cmaxc = whole_props[matched_whole]["bbox"]
            minr = min(nminr, cminr)
            minc = min(nminc, cminc)
            maxr = max(nmaxr, cmaxr)
            maxc = max(nmaxc, cmaxc)
            used_whole.add(matched_whole)

            rs = slice(minr, maxr)
            cs = slice(minc, maxc)
            npatch = nuc[rs, cs] == nlab    # nucleus pixels in the crop
            cpatch = whole[rs, cs] == matched_whole  # cell pixels in the crop

            # Avoid overlaps in synthesized output labels: only paint into
            # pixels not yet claimed by an earlier cell (first-come priority).
            available = out_cell[rs, cs] == 0
            cpatch = cpatch & available
            npatch = npatch & cpatch  # nucleus must also be within available cell pixels

            if not np.any(cpatch):
                dropped_synth_cells += 1
                continue

            # Paint the matched cell and its nucleus into the output arrays
            # with a shared cell_id so they can be cross-referenced later.
            out_cell[rs, cs][cpatch] = next_id
            out_nuc[rs, cs][npatch] = next_id

            # Compute tight bounding box (with 1px padding for contour safety)
            # from the actual painted pixels, not from the input region props.
            cell_rows, cell_cols = np.nonzero(cpatch)
            cell_minr = minr + int(cell_rows.min())
            cell_minc = minc + int(cell_cols.min())
            cell_maxr = minr + int(cell_rows.max()) + 1
            cell_maxc = minc + int(cell_cols.max()) + 1
            bbox_map[next_id] = (
                slice(max(0, cell_minr - 1), min(nuc.shape[0], cell_maxr + 1)),
                slice(max(0, cell_minc - 1), min(nuc.shape[1], cell_maxc + 1)),
            )
            records.append(CellRecord(next_id, matched_whole, int(nlab)))
            next_id += 1
        else:
            # Reserve an id and record it; actual pixels assigned in pass 2.
            unmatched_nuclei.append((nlab, next_id))
            next_id += 1

    # --- Pass 2: watershed-partition unmatched nuclei so they don't overlap ---
    if unmatched_nuclei:
        # Build a seed image: each unmatched nucleus gets its reserved cell id.
        seeds = np.zeros_like(out_cell)
        for nlab, cid in unmatched_nuclei:
            seeds[nuc == nlab] = cid

        # Restrict growth to pixels within estimate_dist of any unmatched nucleus
        # and not already claimed by pass-1 cells.
        any_unmatched_nuc = seeds > 0
        growth_zone = ndi.binary_dilation(
            any_unmatched_nuc,
            structure=disk(max(1, int(round(estimate_dist)))),
        )
        growth_zone = growth_zone & (out_cell == 0)

        # Euclidean distance from each background pixel to the nearest seed;
        # watershed floods lowest-distance pixels first → Voronoi-like partition.
        dist_map = ndi.distance_transform_edt(seeds == 0)
        ws = watershed(dist_map, markers=seeds, mask=growth_zone)

        for nlab, cid in unmatched_nuclei:
            cpatch_full = ws == cid
            npatch_full = (nuc == nlab) & cpatch_full

            if not np.any(cpatch_full):
                dropped_synth_cells += 1
                continue

            out_cell[cpatch_full] = cid
            out_nuc[npatch_full] = cid

            cell_rows, cell_cols = np.nonzero(cpatch_full)
            cell_minr = int(cell_rows.min())
            cell_minc = int(cell_cols.min())
            cell_maxr = int(cell_rows.max()) + 1
            cell_maxc = int(cell_cols.max()) + 1
            bbox_map[cid] = (
                slice(max(0, cell_minr - 1), min(nuc.shape[0], cell_maxr + 1)),
                slice(max(0, cell_minc - 1), min(nuc.shape[1], cell_maxc + 1)),
            )
            records.append(CellRecord(cid, None, int(nlab)))

    unmatched_whole = len(set(whole_labels) - used_whole)
    stats = {
        "nucleus_count": len(nuc_props),
        "whole_cell_count": len(whole_props),
        "matched_cells": len(records),
        "unmatched_whole_cells": unmatched_whole,
        "dropped_synth_cells": dropped_synth_cells,
    }
    return out_cell, out_nuc, records, stats, bbox_map


def mask_to_geometry(
    mask: np.ndarray,
    simplify: bool,
    tolerance: float,
    row_offset: int = 0,
    col_offset: int = 0,
):
    """Convert a binary mask to a Shapely Polygon via marching-squares contours.

    Extracts sub-pixel contours from ``mask`` using ``skimage.find_contours``
    at the 0.5 iso-level (the midpoint between 0 and 1), converts them to
    Shapely Polygons, merges any disjoint pieces via ``unary_union``, and
    optionally simplifies the result with Douglas-Peucker.

    Contour coordinates are shifted by ``(row_offset, col_offset)`` so that
    geometries from local bounding-box crops map back to global image
    coordinates.

    Non-polygon geometry artefacts (LineStrings, Points) that can arise from
    degenerate contours or boolean operations are discarded; only the single
    largest Polygon is returned, matching QuPath's convention of one polygon
    per cell detection.

    Parameters
    ----------
    mask : np.ndarray
        2-D binary mask (truthy = foreground).
    simplify : bool
        Whether to apply Douglas-Peucker simplification.
    tolerance : float
        Simplification tolerance in pixels.  Lower values preserve more shape
        detail at the cost of larger GeoJSON output.
    row_offset, col_offset : int
        Global pixel offsets added to contour coordinates.

    Returns
    -------
    shapely.geometry.Polygon or None
        The cell polygon in global image coordinates, or ``None`` if the
        mask is empty or produces no valid geometry.
    """
    if not np.any(mask):
        return None

    contours = find_contours(mask.astype(np.uint8), level=0.5)
    if not contours:
        return None

    polys = []
    for c in contours:
        if len(c) < 3:
            continue
        xy = [(float(col_offset + p[1]), float(row_offset + p[0])) for p in c]
        poly = Polygon(xy)
        if poly.is_empty:
            continue
        if not poly.is_valid:
            poly = poly.buffer(0)
        if not poly.is_empty:
            polys.append(poly)

    if not polys:
        return None

    g = unary_union(polys)
    if g.geom_type == "GeometryCollection":
        keep = [geom for geom in g.geoms if geom.geom_type in ("Polygon", "MultiPolygon") and not geom.is_empty]
        if not keep:
            return None
        g = unary_union(keep)

    if simplify and tolerance > 0:
        g = g.simplify(tolerance, preserve_topology=True)
    if not g.is_valid:
        g = g.buffer(0)
    if g.is_empty:
        return None

    # Ensure we only return Polygon (not GeometryCollection, MultiPolygon, etc.)
    g = _ensure_largest_polygon(g)
    if g is None or g.is_empty:
        return None

    return g


def basic_shape_metrics(cell_mask: np.ndarray, nuc_mask: Optional[np.ndarray], px_um: float) -> Dict[str, float]:
    """Compute morphological shape descriptors for a single cell and its nucleus.

    Uses ``skimage.measure.regionprops`` on the binary mask to extract:

    * **Area** — pixel count converted to µm² via ``(pixel_size)²``.
    * **Circularity** — ``4π · area / perimeter²``; equals 1.0 for a
      perfect circle, <1 for elongated or irregular shapes.
    * **Length** — perimeter in µm.
    * **Max / Min diameter** — major and minor axis lengths of the
      best-fit ellipse.
    * **Solidity** — ``area / convex_hull_area``; how 'filled' the shape
      is (1.0 = fully convex).

    If a nuclear mask is provided and non-empty, the same metrics are
    computed for the nucleus, plus a **Nucleus/Cell area ratio**.

    .. note::

       These measurements are computed from the **raw segmentation mask**
       before any polygon overlap clipping.  They reflect the segmentation
       model's output, not the final QuPath-compatible geometry.

    Parameters
    ----------
    cell_mask : np.ndarray
        Boolean or binary mask of the cell body.
    nuc_mask : np.ndarray or None
        Boolean or binary mask of the nucleus (same shape as *cell_mask*).
    px_um : float
        Pixel size in micrometres.

    Returns
    -------
    dict
        Measurement name → float value, keyed like
        ``"Cell: Area µm^2"``, ``"Nucleus: Circularity"``, etc.
    """
    rp = regionprops(cell_mask.astype(np.uint8))
    if not rp:
        return {}
    r = rp[0]
    area_um2 = float(r.area) * (px_um ** 2)
    perimeter_um = float(r.perimeter) * px_um if r.perimeter > 0 else 0.0
    circularity = float(4 * math.pi * r.area / (r.perimeter ** 2)) if r.perimeter > 0 else 0.0
    maj = float(r.major_axis_length) * px_um
    mino = float(r.minor_axis_length) * px_um
    solidity = float(r.solidity) if r.solidity is not None else 0.0

    out = {
        "Cell: Area µm^2": area_um2,
        "Cell: Circularity": circularity,
        "Cell: Length µm": perimeter_um,
        "Cell: Max diameter µm": maj,
        "Cell: Min diameter µm": mino,
        "Cell: Solidity": solidity,
    }

    if nuc_mask is not None and np.any(nuc_mask):
        nrp = regionprops(nuc_mask.astype(np.uint8))
        if nrp:
            nr = nrp[0]
            n_area = float(nr.area) * (px_um ** 2)
            out["Nucleus: Area µm^2"] = n_area
            out["Nucleus: Circularity"] = float(4 * math.pi * nr.area / (nr.perimeter ** 2)) if nr.perimeter > 0 else 0.0
            out["Nucleus: Length µm"] = float(nr.perimeter) * px_um if nr.perimeter > 0 else 0.0
            out["Nucleus: Max diameter µm"] = float(nr.major_axis_length) * px_um
            out["Nucleus: Min diameter µm"] = float(nr.minor_axis_length) * px_um
            out["Nucleus: Solidity"] = float(nr.solidity) if nr.solidity is not None else 0.0
            out["Nucleus/Cell area ratio"] = float(n_area / area_um2) if area_um2 > 0 else 0.0

    return out


def compartment_masks(cell_mask: np.ndarray, nuc_mask: np.ndarray) -> Dict[str, np.ndarray]:
    """Derive the four sub-cellular compartment masks from cell and nucleus masks.

    Compartments follow QuPath's cell measurement model:

    * **CELL** — the entire cell body (boolean of *cell_mask*).
    * **NUCLEUS** — the nuclear region (boolean of *nuc_mask*).
    * **CYTOPLASM** — cell minus nucleus: ``cell & ~nucleus``.
    * **MEMBRANE** — 1-pixel-wide outer ring of the cell, obtained by
      subtracting the eroded cell from the original:
      ``cell & ~erode(cell, disk(1))``.  This is a morphological
      approximation of the plasma membrane.

    Parameters
    ----------
    cell_mask : np.ndarray
        Binary cell body mask.
    nuc_mask : np.ndarray
        Binary nucleus mask (same shape).

    Returns
    -------
    dict
        ``{"CELL": ..., "NUCLEUS": ..., "CYTOPLASM": ..., "MEMBRANE": ...}``
        Each value is a boolean ``np.ndarray``.
    """
    cm = cell_mask.astype(bool)
    nm = nuc_mask.astype(bool)
    cyto = cm & ~nm
    mem = cm & ~binary_erosion(cm, _DISK_1)
    return {
        "CELL": cm,
        "NUCLEUS": nm,
        "CYTOPLASM": cyto,
        "MEMBRANE": mem,
    }


def stat_values(vals: np.ndarray) -> Dict[str, float]:
    """Compute summary statistics for a 1-D array of pixel intensity values.

    Returns mean, median, min, max, and standard deviation.  These match
    the summary statistics that QuPath reports per-channel per-compartment.
    Returns an empty dict if *vals* is empty (e.g. the compartment mask
    has zero pixels).
    """
    if vals.size == 0:
        return {}
    return {
        "Mean": float(np.mean(vals)),
        "Median": float(np.median(vals)),
        "Min": float(np.min(vals)),
        "Max": float(np.max(vals)),
        "Std.Dev.": float(np.std(vals)),
    }


def add_intensity_measurements(props: Dict[str, Any], image_cyx: np.ndarray, ch_names: Sequence[str], comp_masks: Dict[str, np.ndarray]):
    """Add per-channel, per-compartment intensity summary statistics to *props*.

    For every combination of image channel and compartment mask (Cell,
    Nucleus, Cytoplasm, Membrane), extracts the pixel values under the
    mask and computes Mean, Median, Min, Max, and Std.Dev.  Keys are
    formatted as ``"<channel>: <compartment>: <stat>"`` to match QuPath's
    measurement table convention.

    Parameters
    ----------
    props : dict
        Measurement dictionary to populate (mutated in place).
    image_cyx : np.ndarray
        Multi-channel image crop, shape ``(C, H, W)``.
    ch_names : sequence of str
        Channel names, length ``C``.
    comp_masks : dict
        Compartment name → boolean mask, as returned by
        ``compartment_masks()``.
    """
    labels = {"CELL": "Cell", "NUCLEUS": "Nucleus", "CYTOPLASM": "Cytoplasm", "MEMBRANE": "Membrane"}
    for ci, ch in enumerate(ch_names):
        ch_img = image_cyx[ci]
        for comp, m in comp_masks.items():
            vals = ch_img[m]
            if vals.size == 0:
                continue
            for k, v in stat_values(vals).items():
                props[f"{ch}: {labels[comp]}: {k}"] = v


def add_percentiles(props: Dict[str, Any], image_cyx: np.ndarray, ch_names: Sequence[str], comp_masks: Dict[str, np.ndarray], percentiles: Sequence[float]):
    """Add user-specified intensity percentiles per channel and compartment.

    Similar to ``add_intensity_measurements`` but computes arbitrary
    percentiles (e.g. 5th, 25th, 75th, 95th) instead of fixed summary
    statistics.  Useful for robust central-tendency estimators or for
    detecting bimodal intensity distributions within a compartment.

    Keys are formatted as ``"<channel>: <compartment>: Percentile: <p>"``.

    Parameters
    ----------
    props : dict
        Measurement dictionary (mutated in place).
    image_cyx : np.ndarray
        Multi-channel image crop, shape ``(C, H, W)``.
    ch_names : sequence of str
        Channel names.
    comp_masks : dict
        Compartment → boolean mask.
    percentiles : sequence of float
        Percentile values to compute (0–100).
    """
    labels = {"CELL": "Cell", "NUCLEUS": "Nucleus", "CYTOPLASM": "Cytoplasm", "MEMBRANE": "Membrane"}
    for ci, ch in enumerate(ch_names):
        ch_img = image_cyx[ci]
        for comp, m in comp_masks.items():
            vals = ch_img[m]
            if vals.size == 0:
                continue
            for p in percentiles:
                props[f"{ch}: {labels[comp]}: Percentile: {p}"] = float(np.percentile(vals, p))


def _erosion_bins_for_mask(mask: np.ndarray, n_bins: int = 5) -> List[Tuple[np.ndarray, int]]:
    """Compute erosion depths that divide a binary mask into n equal-area bins.

    Works inward from the cell boundary. Each bin targets 1/n_bins of the
    total mask area. Returns a list of (eroded_mask, depth_px) tuples, one
    per bin, where each mask represents the region at that depth.

    The approach: iteratively erode by 1 pixel at a time, recording the
    cumulative area removed. When cumulative removal crosses the next
    (bin_index / n_bins) * total_area threshold, we snapshot that erosion
    level as the bin boundary.
    """
    total = int(np.count_nonzero(mask))
    if total == 0:
        return []

    target_fractions = [(b / n_bins) for b in range(1, n_bins + 1)]  # 0.2, 0.4, 0.6, 0.8, 1.0
    bins: List[Tuple[np.ndarray, int]] = []

    current = mask.astype(bool)
    depth = 0

    for target_frac in target_fractions:
        target_remaining = int(total * (1.0 - target_frac))  # area we want left after this bin boundary
        # Erode until area drops to or below target_remaining, or mask empties
        while True:
            area = int(np.count_nonzero(current))
            if area <= target_remaining or area == 0:
                break
            current = binary_erosion(current, _DISK_1)
            depth += 1
        bins.append((current.copy(), depth))
        if area == 0:
            # All remaining bins will also be empty — pad with empty masks
            while len(bins) < n_bins:
                bins.append((current.copy(), depth))
            break

    return bins


def add_erosion_measurements(
    props: Dict[str, Any],
    image_cyx: np.ndarray,
    ch_names: Sequence[str],
    comp_masks: Dict[str, np.ndarray],
    steps: Sequence[int],  # kept for API compatibility but ignored
    n_bins: int = 5,
):
    """Measure intensity in 5 concentric erosion bins, each covering ~20% of cell area.

    Instead of fixed pixel depths, bins are defined by equal area fractions
    working inward from the cell boundary. Bin 1 is the outermost 20% ring,
    Bin 5 is the innermost 20% (deepest core). This makes bins comparable
    across cells of different sizes.

    The ``steps`` parameter is accepted for API compatibility but ignored;
    bin boundaries are computed adaptively from the mask geometry.

    Produces per bin:
      * ``<Compartment>: ErosionBin_<N>: Area_px`` — absolute pixel count
      * ``<Compartment>: ErosionBin_<N>: Area_Fraction`` — fraction of total compartment area
      * ``<Compartment>: ErosionBin_<N>: Depth_px`` — erosion depth at the inner edge of this bin
      * ``<channel>: <Compartment>: ErosionBin_<N>: Mean``
      * ``<channel>: <Compartment>: ErosionBin_<N>: Median``
    """
    for comp in ("CELL", "NUCLEUS"):
        base = comp_masks[comp]
        base_area = int(np.count_nonzero(base))
        if base_area == 0:
            continue

        comp_name = comp.capitalize()
        bin_boundaries = _erosion_bins_for_mask(base, n_bins=n_bins)

        # Convert cumulative eroded masks → mutually exclusive annular rings
        # Ring N = (area remaining after bin N-1) minus (area remaining after bin N)
        prev_mask = base.astype(bool)
        for bin_idx, (eroded_mask, depth_px) in enumerate(bin_boundaries, start=1):
            ring = prev_mask & ~eroded_mask  # the shell peeled off in this bin
            ring_area = int(np.count_nonzero(ring))

            props[f"{comp_name}: ErosionBin_{bin_idx}: Area_px"] = ring_area
            props[f"{comp_name}: ErosionBin_{bin_idx}: Area_Fraction"] = float(ring_area / base_area)
            props[f"{comp_name}: ErosionBin_{bin_idx}: Depth_px"] = depth_px

            if ring_area > 0:
                for ci, ch in enumerate(ch_names):
                    vals = image_cyx[ci][ring]
                    if vals.size > 0:
                        props[f"{ch}: {comp_name}: ErosionBin_{bin_idx}: Mean"] = float(np.mean(vals))
                        props[f"{ch}: {comp_name}: ErosionBin_{bin_idx}: Median"] = float(np.median(vals))

            prev_mask = eroded_mask


def _expansion_bins_for_mask(cell_mask: np.ndarray, total_expansion_px: int, n_bins: int = 5) -> List[Tuple[np.ndarray, int]]:
    """Compute dilation depths that divide a 20 µm expansion zone into n equal-area bins.

    Works outward from the cell boundary. First dilates the cell mask by
    *total_expansion_px* pixels to define the full expansion zone, then
    iteratively dilates from the cell boundary 1 pixel at a time, snapshotting
    when cumulative ring area crosses each (bin_index / n_bins) * total_zone_area
    threshold.

    Returns a list of (dilated_mask, depth_px) tuples, one per bin.
    """
    cm = cell_mask.astype(bool)
    if not np.any(cm):
        return []

    full_dilated = ndi.binary_dilation(cm, structure=_DISK_1, iterations=total_expansion_px)
    zone = full_dilated & ~cm
    total_zone_area = int(np.count_nonzero(zone))
    if total_zone_area == 0:
        return []

    target_fractions = [(b / n_bins) for b in range(1, n_bins + 1)]  # 0.2, 0.4, 0.6, 0.8, 1.0
    bins: List[Tuple[np.ndarray, int]] = []

    current = cm.copy()
    depth = 0

    for target_frac in target_fractions:
        target_area = int(total_zone_area * target_frac)
        # Dilate until cumulative ring area reaches target
        while depth < total_expansion_px:
            current_ring_area = int(np.count_nonzero(current & ~cm))
            if current_ring_area >= target_area:
                break
            current = ndi.binary_dilation(current, structure=_DISK_1, iterations=1)
            depth += 1
        bins.append((current.copy(), depth))
        if depth >= total_expansion_px:
            while len(bins) < n_bins:
                bins.append((current.copy(), depth))
            break

    return bins


def add_expansion_measurements(
    props: Dict[str, Any],
    image_cyx: np.ndarray,
    ch_names: Sequence[str],
    cell_mask: np.ndarray,
    pixel_size_microns: float,
    n_bins: int = 5,
):
    """Measure intensity in 5 equal-area annular bins dilated 20 µm outward from the cell.

    Instead of fixed pixel steps, the full 20 µm expansion zone is divided
    into *n_bins* rings of approximately equal area. Bin 1 is the ring
    immediately adjacent to the cell boundary, Bin 5 is the outermost ring
    at the edge of the 20 µm zone. This makes bins comparable across cells
    of different sizes and across imaging platforms with different pixel sizes.

    Produces per bin:
      * ``Cell: ExpansionBin_<N>: Area_px`` — absolute pixel count
      * ``Cell: ExpansionBin_<N>: Area_Fraction`` — fraction of cell body area
      * ``Cell: ExpansionBin_<N>: Depth_px`` — dilation depth at the outer edge
      * ``<channel>: Cell: ExpansionBin_<N>: Mean``
      * ``<channel>: Cell: ExpansionBin_<N>: Median``

    Parameters
    ----------
    props : dict
        Measurement dictionary (mutated in place).
    image_cyx : np.ndarray
        Multi-channel image crop.
    ch_names : sequence of str
        Channel names.
    cell_mask : np.ndarray
        Binary cell body mask.
    pixel_size_microns : float
        Pixel size in microns, used to convert 20 µm to pixels.
    n_bins : int
        Number of equal-area bins (default 5).
    """
    EXPANSION_UM = 20.0
    total_expansion_px = max(1, int(round(EXPANSION_UM / pixel_size_microns)))

    cm = cell_mask.astype(bool)
    base_area = int(np.count_nonzero(cm))
    if base_area == 0:
        return

    bin_boundaries = _expansion_bins_for_mask(cm, total_expansion_px, n_bins=n_bins)
    if not bin_boundaries:
        return

    # Convert cumulative dilated masks → mutually exclusive annular rings
    prev_mask = cm.copy()
    for bin_idx, (dilated_mask, depth_px) in enumerate(bin_boundaries, start=1):
        ring = dilated_mask & ~prev_mask
        ring_area = int(np.count_nonzero(ring))

        props[f"Cell: ExpansionBin_{bin_idx}: Area_px"] = ring_area
        props[f"Cell: ExpansionBin_{bin_idx}: Area_Fraction"] = float(ring_area / base_area)
        props[f"Cell: ExpansionBin_{bin_idx}: Depth_px"] = depth_px

        if ring_area > 0:
            for ci, ch in enumerate(ch_names):
                vals = image_cyx[ci][ring]
                if vals.size > 0:
                    props[f"{ch}: Cell: ExpansionBin_{bin_idx}: Mean"] = float(np.mean(vals))
                    props[f"{ch}: Cell: ExpansionBin_{bin_idx}: Median"] = float(np.median(vals))

        prev_mask = dilated_mask


def add_environment_measurements(
    props: Dict[str, Any],
    image_cyx: np.ndarray,
    ch_names: Sequence[str],
    cell_mask: np.ndarray,
    pixel_size_microns: float,
):
    """Measure intensity in a 20 µm pericellular environment zone.

    Unlike ``add_expansion_measurements`` which splits the dilated zone into
    mutually exclusive annular rings, this function computes a single
    dilation of 20 µm (converted to pixels via *pixel_size_microns*) and
    measures the **entire** zone between the cell boundary and the outer
    ring as one compartment.  This matches CellTune's "Environment"
    compartment for characterising the local tissue microenvironment
    around each cell.

    Produces:
      * ``Cell: Environment_20um: Pixel_Count`` and ``Area_Fraction``
      * ``<channel>: Cell: Environment_20um: Mean/Median/Min/Max/Std.Dev.``

    Parameters
    ----------
    props : dict
        Measurement dictionary (mutated in place).
    image_cyx : np.ndarray
        Multi-channel image crop.
    ch_names : sequence of str
        Channel names.
    cell_mask : np.ndarray
        Binary cell body mask.
    pixel_size_microns : float
        Pixel size in microns, used to convert the 20 µm radius to pixels.
    """
    ENVIRONMENT_UM = 20.0
    expansion_px = max(1, int(round(ENVIRONMENT_UM / pixel_size_microns)))
    cm = cell_mask.astype(bool)
    if not np.any(cm):
        return
    dilated = ndi.binary_dilation(cm, structure=_DISK_1, iterations=expansion_px)
    env_mask = dilated & ~cm
    env_area = int(np.count_nonzero(env_mask))
    if env_area == 0:
        return
    base_area = int(np.count_nonzero(cm))
    props["Cell: Environment_20um: Pixel_Count"] = env_area
    props["Cell: Environment_20um: Area_Fraction"] = float(env_area / base_area) if base_area > 0 else 0.0
    for ci, ch in enumerate(ch_names):
        vals = image_cyx[ci][env_mask]
        if vals.size == 0:
            continue
        props[f"{ch}: Cell: Environment_20um: Mean"] = float(np.mean(vals))
        props[f"{ch}: Cell: Environment_20um: Median"] = float(np.median(vals))
        props[f"{ch}: Cell: Environment_20um: Min"] = float(np.min(vals))
        props[f"{ch}: Cell: Environment_20um: Max"] = float(np.max(vals))
        props[f"{ch}: Cell: Environment_20um: Std.Dev."] = float(np.std(vals))


def add_neighborhood_features(features: List[Dict[str, Any]], k: int, pixel_size_microns: float = 0.5) -> None:
    """Aggregate each cell's numeric measurements across its k nearest neighbours
    within a 20 µm radius. Only the mean is computed (max is omitted).

    Cells with no neighbours within the distance cap produce no neighbourhood
    keys, so isolated cells in sparse tissue don't get meaningless aggregations
    from cells hundreds of microns away.
    """
    #TODO, is mean enough? Worth adding median and other statistics?
    if k <= 0 or len(features) < 2:
        return

    MAX_DISTANCE_UM = 20.0
    max_distance_px = MAX_DISTANCE_UM / pixel_size_microns

    centroids = np.empty((len(features), 2), dtype=np.float64)
    for i, feat in enumerate(features):
        geom = shape(feat["geometry"])
        c = geom.centroid
        centroids[i] = (c.x, c.y)

    tree = cKDTree(centroids)
    actual_k = min(k + 1, len(features))
    distances, indices = tree.query(centroids, k=actual_k)

    sample_meas = features[0]["properties"].get("measurements", {})
    numeric_keys = [
        key for key, val in sample_meas.items()
        if isinstance(val, (int, float)) and key not in ("id", "cell_label", "nucleus_label")
    ]
    if not numeric_keys:
        return

    n = len(features)
    print(f"  Extracting {len(numeric_keys)} measurement vectors for {n} cells (20 µm cap = {max_distance_px:.1f} px)...")
    meas_vectors: Dict[str, np.ndarray] = {}
    for ki, key in enumerate(numeric_keys):
        arr = np.full(n, np.nan, dtype=np.float64)
        for i, feat in enumerate(features):
            v = feat["properties"].get("measurements", {}).get(key)
            if v is not None:
                arr[i] = v
        meas_vectors[key] = arr
        if (ki + 1) % 100 == 0 or (ki + 1) == len(numeric_keys):
            print(f"  Extracted {ki + 1}/{len(numeric_keys)} measurement vectors")

    for i, feat in enumerate(features):
        if (i + 1) % 10000 == 0 or (i + 1) == n:
            print(f"  Neighborhood aggregation: {i + 1}/{n} cells ({int((i + 1) * 100.0 / max(n, 1))}%)")
        if actual_k <= 1:
            continue

        nbr_idx = indices[i, 1:actual_k]
        nbr_dist = distances[i, 1:actual_k]

        # Apply 20 µm distance cap — drop neighbours beyond the threshold
        within = nbr_dist <= max_distance_px
        nbr_idx = nbr_idx[within]

        if len(nbr_idx) == 0:
            continue  # isolated cell — skip rather than writing NaN features

        meas = feat["properties"].setdefault("measurements", {})
        for key in numeric_keys:
            vals = meas_vectors[key][nbr_idx]
            valid = vals[~np.isnan(vals)]
            if valid.size == 0:
                continue
            meas[f"Neighbors: Mean: {key}"] = float(np.mean(valid))


def parse_csv_numbers(s: str, cast=float, positive_only=False) -> List:
    """Parse a comma-separated string of numbers into a typed list.

    Used to parse ``--percentiles`` CLI arguments.
    Empty or whitespace-only strings return an empty list.  Invalid tokens
    raise ``ValueError``.

    Parameters
    ----------
    s : str
        Comma-separated numeric string, e.g. ``"5,25,75,95"``.
    cast : callable
        Type constructor (``float`` or ``int``).
    positive_only : bool
        If True, silently drop values ≤ 0.

    Returns
    -------
    list
        Parsed and optionally filtered numeric values.
    """
    if not s or not s.strip():
        return []
    out = []
    for x in s.split(","):
        x = x.strip()
        if not x:
            continue
        try:
            v = cast(x)
        except (TypeError, ValueError) as e:
            raise ValueError(f"Invalid value '{x}' in list '{s}'") from e
        if positive_only and v <= 0:
            continue
        out.append(v)
    return out


def feature_for_cell(
    cell_id: int,
    cell_crop: np.ndarray,
    nuc_crop: np.ndarray,
    image_crop: np.ndarray,
    row_offset: int,
    col_offset: int,
    rec_cell_label: Optional[int],
    rec_nucleus_label: Optional[int],
    simplify_rois: bool,
    tolerance: float,
    pixel_size_microns: float,
    skip_measurements: bool,
    ch_names: Sequence[str],
    percentiles: Sequence[float],
    erosion_enabled: bool = False,
    expansion_enabled: bool = False,
    environment_expansion: bool = False,
):
    """Build a complete GeoJSON Feature dict for a single cell.

    This is the per-cell workhorse called once per cell (possibly in a
    worker process when ``--threads > 1``).  It receives a small crop of
    the label masks and image around the cell's bounding box and:

    1. Extracts the binary cell mask (``cell_crop == cell_id``) and
       nuclear mask from the label crops.
    2. Converts the cell mask to a Shapely Polygon via ``mask_to_geometry``.
    3. Computes all requested measurements from the **raw raster mask**:
       shape metrics, per-channel per-compartment intensity statistics,
       optional percentiles, erosion profiles, expansion rings, and
       environment zone.
    4. Packages everything into a GeoJSON Feature with ``geometry``,
       ``nucleusGeometry``, and ``properties.measurements``.

    .. note::

       Measurements are computed from the original segmentation mask,
       *before* any polygon overlap clipping in ``constrain_cell_overlaps``.
       This means the measurement values reflect the segmentation model's
       raw output rather than the final display geometry.

    Parameters
    ----------
    cell_id : int
        Unified cell ID (label value in the output mask arrays).
    cell_crop, nuc_crop : np.ndarray
        Label mask crops around this cell's bounding box.
    image_crop : np.ndarray
        Multi-channel image crop, shape ``(C, h, w)``.
    row_offset, col_offset : int
        Global offsets to translate crop coordinates back to full-image
        coordinates for the output geometry.
    rec_cell_label, rec_nucleus_label : int or None
        Original label values from the input masks (for provenance).
    simplify_rois : bool
        Whether to apply Douglas-Peucker simplification to the polygon.
    tolerance : float
        Simplification tolerance in pixels.
    pixel_size_microns : float
        Pixel size for converting area/length to physical units.
    skip_measurements : bool
        If True, only shape metrics are computed (no intensity stats).
    ch_names : sequence of str
        Channel names.
    percentiles, erosion_enabled : sequences/bool
        Optional measurement parameters.
    expansion_enabled : bool
        Whether to compute 20 µm expansion bin measurements.
    environment_expansion : bool
        Whether to compute 20 µm pericellular environment measurements.

    Returns
    -------
    dict or None
        GeoJSON Feature dict, or ``None`` if the cell mask produces no
        valid geometry (e.g. entirely occluded by prior cells).
    """
    # --- Step 1: Extract binary masks from the label crop ---
    cmask = cell_crop == cell_id
    nmask = nuc_crop == cell_id

    # --- Step 2: Convert raster mask to vector polygon (global coords) ---
    geom = mask_to_geometry(cmask, simplify_rois, tolerance, row_offset=row_offset, col_offset=col_offset)
    if geom is None:
        return None
    nuc_geom = mask_to_geometry(nmask, simplify_rois, tolerance, row_offset=row_offset, col_offset=col_offset)

    # --- Step 3: Compute measurements from the RAW raster mask ---
    # Note: these are intentionally computed before overlap clipping so they
    # reflect the segmentation model's original output, not the QuPath-
    # compatible trimmed polygons.
    measurements: Dict[str, Any] = {}
    measurements.update(basic_shape_metrics(cmask, nmask, pixel_size_microns))

    if not skip_measurements:
        # Derive sub-cellular compartment masks (cell, nucleus, cytoplasm, membrane)
        comps = compartment_masks(cmask, nmask)
        add_intensity_measurements(measurements, image_crop, ch_names, comps)
        if percentiles:
            add_percentiles(measurements, image_crop, ch_names, comps, percentiles)
        if erosion_enabled:
            add_erosion_measurements(measurements, image_crop, ch_names, comps, steps=[])
        if expansion_enabled:
            add_expansion_measurements(measurements, image_crop, ch_names, cmask, pixel_size_microns)
        if environment_expansion:
            add_environment_measurements(measurements, image_crop, ch_names, cmask, pixel_size_microns)

    # --- Step 4: Package into GeoJSON Feature ---
    feature: Dict[str, Any] = {
        "type": "Feature",
        "id": f"cell-{cell_id}",
        "geometry": mapping(geom),
        "properties": {
            "objectType": "cell",
            "id": int(cell_id),
            "cell_label": int(rec_cell_label) if rec_cell_label is not None else None,
            "nucleus_label": int(rec_nucleus_label) if rec_nucleus_label is not None else None,
            "measurements": measurements,
        },
    }
    if nuc_geom is not None:
        feature["nucleusGeometry"] = mapping(nuc_geom)

    return feature


def iter_tasks(
    unique_cells: Sequence[int],
    bbox_map: Dict[int, Tuple[slice, slice]],
    records_by_id: Dict[int, CellRecord],
    cell_labels: np.ndarray,
    nuc_labels: np.ndarray,
    img_cyx: np.ndarray,
    args: argparse.Namespace,
    ch_names: Sequence[str],
    percentiles: Sequence[float],
    erosion_enabled: bool = False,
    expansion_enabled: bool = False,
    environment_expansion: bool = False,
):
    """Lazily yield per-cell task tuples for ``feature_for_cell``.

    Each yielded tuple contains all the data a worker process needs to
    compute one cell's feature independently: a crop of the label masks
    and image around the cell's bounding box, plus all scalar parameters.

    Crops are **copied** from the parent arrays before yielding so that
    each task is self-contained and safely picklable for
    ``ProcessPoolExecutor``.  The generator is consumed lazily to avoid
    materialising all crops in memory simultaneously.

    The bounding box is padded outward by the 20 µm expansion distance
    (in pixels) when expansion or environment measurements are requested,
    so that dilation in the worker doesn't run past the crop boundary.

    Parameters
    ----------
    unique_cells : sequence of int
        Sorted list of cell IDs to process.
    bbox_map : dict
        ``{cell_id: (row_slice, col_slice)}`` from ``match_cells``.
    records_by_id : dict
        ``{cell_id: CellRecord}`` for provenance labels.
    cell_labels, nuc_labels : np.ndarray
        Full-image label masks.
    img_cyx : np.ndarray
        Full multi-channel image.
    args : argparse.Namespace
        Parsed CLI arguments.
    ch_names, percentiles
        Measurement configuration.
    erosion_enabled, expansion_enabled, environment_expansion : bool
        Whether to compute erosion / expansion / environment zone measurements.

    Yields
    ------
    tuple
        Positional arguments for ``feature_for_cell(*)``.
    """
    expand_px = max(1, int(round(20.0 / args.pixel_size_microns))) if (expansion_enabled or environment_expansion) else 0
    max_expand = expand_px
    h, w = cell_labels.shape[:2]
    for cid in unique_cells:
        rs, cs = bbox_map[cid]
        if max_expand > 0:
            r0 = max(rs.start - max_expand, 0)
            r1 = min(rs.stop + max_expand, h)
            c0 = max(cs.start - max_expand, 0)
            c1 = min(cs.stop + max_expand, w)
            rs_pad = slice(r0, r1)
            cs_pad = slice(c0, c1)
        else:
            rs_pad, cs_pad = rs, cs
        rec = records_by_id.get(cid)
        # Copy in parent process to keep each task payload self-contained for worker pickling.
        yield (
            cid,
            cell_labels[rs_pad, cs_pad].copy(),
            nuc_labels[rs_pad, cs_pad].copy(),
            img_cyx[:, rs_pad, cs_pad].copy(),
            rs_pad.start,
            cs_pad.start,
            rec.cell_label if rec else None,
            rec.nucleus_label if rec else None,
            args.simplify_rois,
            args.tolerance,
            args.pixel_size_microns,
            args.skip_measurements,
            tuple(ch_names),
            tuple(percentiles),
            erosion_enabled,
            expansion_enabled,
            environment_expansion,
        )


def iter_tasks_coords(
    unique_cells: Sequence[int],
    bbox_map: Dict[int, Tuple[slice, slice]],
    records_by_id: Dict[int, CellRecord],
    img_h: int,
    img_w: int,
    args: argparse.Namespace,
    ch_names: Sequence[str],
    percentiles: Sequence[float],
    erosion_enabled: bool = False,
    expansion_enabled: bool = False,
    environment_expansion: bool = False,
):
    """Yield per-cell bounding-box coordinate tuples for ``_feature_for_cell_global``.

    Unlike ``iter_tasks``, **no array data is copied or yielded**.  Only scalar
    integers and booleans are yielded, so pickling task payloads to worker
    processes is essentially free.  Workers read crops directly from the
    fork-inherited module globals ``_GLOBAL_IMG``, ``_GLOBAL_CELL``, and
    ``_GLOBAL_NUC`` via copy-on-write shared memory.
    """
    expand_px = max(1, int(round(20.0 / args.pixel_size_microns))) if (expansion_enabled or environment_expansion) else 0
    for cid in unique_cells:
        rs, cs = bbox_map[cid]
        if expand_px > 0:
            r0 = max(rs.start - expand_px, 0)
            r1 = min(rs.stop + expand_px, img_h)
            c0 = max(cs.start - expand_px, 0)
            c1 = min(cs.stop + expand_px, img_w)
        else:
            r0, r1 = rs.start, rs.stop
            c0, c1 = cs.start, cs.stop
        rec = records_by_id.get(cid)
        yield (
            cid, r0, r1, c0, c1,
            rec.cell_label if rec else None,
            rec.nucleus_label if rec else None,
            args.simplify_rois,
            args.tolerance,
            args.pixel_size_microns,
            args.skip_measurements,
            tuple(ch_names),
            tuple(percentiles),
            erosion_enabled,
            expansion_enabled,
            environment_expansion,
        )


def _feature_for_cell_global(
    cell_id: int,
    r0: int, r1: int, c0: int, c1: int,
    rec_cell_label: Optional[int],
    rec_nucleus_label: Optional[int],
    simplify_rois: bool,
    tolerance: float,
    pixel_size_microns: float,
    skip_measurements: bool,
    ch_names: Sequence[str],
    percentiles: Sequence[float],
    erosion_enabled: bool = False,
    expansion_enabled: bool = False,
    environment_expansion: bool = False,
):
    """Worker entry point that reads crops from fork-inherited global arrays.

    Called by worker processes spawned by ``ProcessPoolExecutor``.  Because
    Linux uses ``fork()`` to create workers, the parent's ``_GLOBAL_IMG``,
    ``_GLOBAL_CELL``, and ``_GLOBAL_NUC`` arrays are accessible in each
    worker via copy-on-write without any data transfer.  Workers only READ
    from these globals so pages are never CoW-copied.

    The image crop is cast to float32 here (not at load time) so the parent
    can hold the image in its original dtype (typically uint16), halving
    memory vs. a process-wide float32 copy (~35 GB saved for a 34-channel
    whole-slide image).
    """
    cell_crop = _GLOBAL_CELL[r0:r1, c0:c1].copy()
    if _GLOBAL_SKIP_NUC:
        nuc_crop = np.zeros((r1 - r0, c1 - c0), dtype=np.int64)
    else:
        nuc_crop = _GLOBAL_NUC[r0:r1, c0:c1].copy()
    # Cast only this small crop to float32 — not the full image.
    img_crop = _GLOBAL_IMG[:, r0:r1, c0:c1].astype(np.float32)
    return feature_for_cell(
        cell_id, cell_crop, nuc_crop, img_crop, r0, c0,
        rec_cell_label, rec_nucleus_label, simplify_rois, tolerance,
        pixel_size_microns, skip_measurements, ch_names, percentiles,
        erosion_enabled, expansion_enabled, environment_expansion,
    )


def _ensure_largest_polygon(geom):
    """Extract the largest Polygon from any geometry, discarding non-polygon parts.

    After boolean operations like ``.difference()`` or ``.intersection()``,
    a Shapely geometry may degenerate into a *GeometryCollection*
    containing a mix of Polygons, LineStrings, and Points.  QuPath's
    GeoJSON import requires pure Polygon or MultiPolygon geometry, so
    this function cleans up the result by:

    * Returning *geom* unchanged if it is already a Polygon.
    * Extracting the **largest** Polygon (by area) from a MultiPolygon
      or GeometryCollection.
    * Returning an empty Polygon for any geometry type that is not
      area-bearing (LineString, Point, etc.).
    * Repairing invalid geometries via ``buffer(0)``.

    This mirrors the Cellpose QuPath extension's
    ``GeometryTools.ensurePolygonal()`` + keep-largest logic.
    """
    if geom is None or geom.is_empty:
        return geom
    if not geom.is_valid:
        geom = geom.buffer(0)
    if geom.geom_type == "Polygon":
        return geom
    if geom.geom_type == "MultiPolygon":
        # Keep only the largest polygon
        largest = max(geom.geoms, key=lambda g: g.area)
        return largest
    if geom.geom_type == "GeometryCollection":
        polys = [g for g in geom.geoms if g.geom_type == "Polygon" and not g.is_empty and g.area > 0]
        if not polys:
            return Polygon()  # empty
        return max(polys, key=lambda g: g.area)
    # LineString, Point, etc. — not usable
    return Polygon()


def constrain_cell_overlaps(features: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Clip overlapping cell geometries so no two cells share area.

    Equivalent to QuPath's ``CellTools.constrainCellOverlaps()``.  When two
    cell polygons overlap, the **larger** cell is trimmed via
    ``Polygon.difference()`` so the overlap region is assigned to the
    **smaller** cell.  This heuristic preserves small cells (which are
    more likely to be fully enclosed) and matches QuPath's behaviour.

    **Algorithm:**

    1. **Broad phase — grid spatial hash.**  Each polygon's bounding box is
       inserted into a uniform grid whose cell size targets ~√n grid cells
       (so each bucket has O(1) polygons on average).  Only pairs sharing
       at least one grid bucket are tested for intersection.
    2. **Narrow phase — pairwise Shapely intersection test.**  For each
       candidate pair, ``intersects()`` and ``intersection()`` determine
       whether they genuinely overlap (area > 1e-10 px²).  If so, the
       larger cell's geometry is replaced with
       ``larger.difference(smaller)``, and cleaned up via
       ``_ensure_largest_polygon()``.
    3. **Output assembly.**  Cells whose geometry was reduced to empty by
       successive clipping are dropped.  Nucleus geometries that now
       extend beyond the trimmed cell polygon are intersected down.

    The ``checked`` set ensures each pair is tested at most once even if
    both polygons span multiple grid buckets.

    .. note::

       This function modifies **geometry only**; the ``measurements`` dict
       on each feature is *not* recalculated.  Measurements were computed
       from the raw raster mask in ``feature_for_cell`` and reflect the
       segmentation model's output, not the clipped display polygon.

    Parameters
    ----------
    features : list of dict
        GeoJSON Feature dicts with ``geometry`` and optionally
        ``nucleusGeometry``.

    Returns
    -------
    list of dict
        Features with non-overlapping cell geometries.  Empty cells are
        removed; the list may be shorter than the input.
    """
    if not features:
        return features

    n = len(features)
    # Deserialize GeoJSON geometries into Shapely objects for spatial ops
    geoms = [shape(f["geometry"]) for f in features]
    areas = [g.area for g in geoms]

    # -- broad-phase: grid spatial hash based on geometry bounding boxes --
    # Instead of testing all n*(n-1)/2 pairs, we partition space into a
    # uniform grid and only test pairs that share at least one grid cell.
    bounds = [g.bounds for g in geoms]  # (minx, miny, maxx, maxy)
    all_minx = min(b[0] for b in bounds)
    all_miny = min(b[1] for b in bounds)
    all_maxx = max(b[2] for b in bounds)
    all_maxy = max(b[3] for b in bounds)
    span = max(all_maxx - all_minx, all_maxy - all_miny, 1.0)
    # Grid cell size is chosen so there are ~sqrt(n) buckets along each
    # axis, meaning each bucket contains ~1 geometry on average for a
    # uniform distribution.  Dense clusters will have more, but the
    # 'checked' set prevents redundant pair tests.
    grid_size = max(span / max(int(n ** 0.5), 1), 1.0)

    print(f"  Building spatial grid (grid_size={grid_size:.1f})...")
    # Insert each geometry into every grid bucket its bounding box overlaps
    grid: Dict[Tuple[int, int], List[int]] = {}
    for i, (minx, miny, maxx, maxy) in enumerate(bounds):
        gx0 = int((minx - all_minx) / grid_size)
        gy0 = int((miny - all_miny) / grid_size)
        gx1 = int((maxx - all_minx) / grid_size)
        gy1 = int((maxy - all_miny) / grid_size)
        for gx in range(gx0, gx1 + 1):
            for gy in range(gy0, gy1 + 1):
                grid.setdefault((gx, gy), []).append(i)

    # -- narrow-phase: pairwise intersection test and clipping --
    # For each grid bucket, test all pairs of geometries within it.
    # The 'checked' set avoids re-testing pairs that appear in multiple buckets.
    checked: set = set()
    clipped = 0
    _grid_total = len(grid)
    _grid_done = 0
    for cell_list in grid.values():
        _grid_done += 1
        if _grid_done % 5000 == 0 or _grid_done == _grid_total:
            print(f"  Overlap constraint progress: {_grid_done}/{_grid_total} grid cells, "
                  f"{len(checked)} pairs checked, {clipped} clipped")
        for ii in range(len(cell_list)):
            i = cell_list[ii]
            for jj in range(ii + 1, len(cell_list)):
                j = cell_list[jj]
                pair = (min(i, j), max(i, j))
                if pair in checked:
                    continue
                checked.add(pair)

                gi = geoms[i]
                gj = geoms[j]
                if gi.is_empty or gj.is_empty:
                    continue

                try:
                    if not gi.intersects(gj):
                        continue
                    intersection = gi.intersection(gj)
                    if intersection.is_empty or intersection.area < 1e-10:
                        continue
                except Exception:
                    continue

                # Trim the larger cell; keep the smaller cell intact.
                if areas[i] >= areas[j]:
                    gi = _ensure_largest_polygon(gi.difference(gj))
                    geoms[i] = gi
                    areas[i] = gi.area if gi is not None and not gi.is_empty else 0
                else:
                    gj = _ensure_largest_polygon(gj.difference(gi))
                    geoms[j] = gj
                    areas[j] = gj.area if gj is not None and not gj.is_empty else 0
                clipped += 1

    print(f"  Narrow phase complete: {len(checked)} pairs checked, {clipped} clipped. Building output...")
    out = []
    for f, g in zip(features, geoms):
        if g is None or g.is_empty:
            continue
        # Final safety: ensure only Polygon/MultiPolygon in output
        g = _ensure_largest_polygon(g)
        if g is None or g.is_empty:
            continue
        f["geometry"] = mapping(g)
        # Also clip nucleusGeometry if it extends beyond the trimmed cell
        if "nucleusGeometry" in f:
            try:
                ng = shape(f["nucleusGeometry"])
                ng = ng.intersection(g)
                ng = _ensure_largest_polygon(ng)
                if ng is not None and not ng.is_empty:
                    f["nucleusGeometry"] = mapping(ng)
                else:
                    del f["nucleusGeometry"]
            except Exception:
                del f["nucleusGeometry"]
        out.append(f)

    print(f"Overlap constraint: checked {len(checked)} pairs, clipped {clipped}, removed {n - len(out)} empty cells")
    return out


def rasterize_features_to_mask(
    features: List[Dict[str, Any]], height: int, width: int
) -> np.ndarray:
    """Rasterize cell feature polygons back to an integer label mask TIFF.

    Each cell is rasterized with its ``properties.id`` as the label value.
    Produces the same format as smooth_masks output — a 2-D integer label
    image where 0 is background.

    This is the inverse of the contour extraction in ``mask_to_geometry``:
    it converts the final post-clipping vector polygons back to a raster
    representation.  The resulting mask reflects the **clipped** geometries
    (after ``constrain_cell_overlaps``), so it is suitable for use in
    downstream raster-based analyses that need non-overlapping cell regions.

    Rasterization uses ``skimage.draw.polygon`` (scan-line fill).  For
    MultiPolygon cells, all component polygons are filled.  Interior holes
    in polygons are **not** handled (they are filled), which matches the
    typical cell segmentation use case where cells are simply-connected.

    Parameters
    ----------
    features : list of dict
        GeoJSON Feature dicts (only ``objectType == "cell"`` are used).
    height, width : int
        Output mask dimensions (should match the original image).

    Returns
    -------
    np.ndarray
        2-D int32 or int64 label mask of shape ``(height, width)``.
    """
    from skimage.draw import polygon as draw_polygon

    max_id = max(
        (f["properties"].get("id", 0) for f in features if f["properties"].get("objectType") == "cell"),
        default=0,
    )
    dtype = np.int32 if max_id < 2**31 else np.int64
    mask = np.zeros((height, width), dtype=dtype)

    for feat in features:
        if feat["properties"].get("objectType") != "cell":
            continue
        cell_id = feat["properties"].get("id", 0)
        if cell_id <= 0:
            continue

        geom = shape(feat["geometry"])
        if geom.is_empty:
            continue

        polys = []
        if geom.geom_type == "Polygon":
            polys = [geom]
        elif geom.geom_type == "MultiPolygon":
            polys = list(geom.geoms)

        for poly in polys:
            ext = np.array(poly.exterior.coords)
            rr, cc = draw_polygon(ext[:, 1], ext[:, 0], shape=(height, width))
            mask[rr, cc] = cell_id

    return mask


def regularize_mask(
    mask: np.ndarray,
    seed_centroids: Optional[Dict[int, Tuple[float, float]]] = None,
) -> np.ndarray:
    """Re-partition a label mask using watershed from seed points.

    When overlapping instance masks (e.g. from CellSAM) are flattened to a
    single label plane, the last-written label 'wins' overlap pixels, creating
    irregular, jagged boundaries between adjacent cells.  This function
    redistributes those pixels using watershed flooding from seed centroids,
    producing smooth, equidistant (Voronoi-like) boundaries.

    **How it works:**

    1. Compute the Euclidean distance transform from each seed point.
    2. Run watershed on this distance image, restricted to the occupied
       (non-background) region of the original mask.
    3. Each pixel is assigned to the label of the nearest seed, producing
       equidistant boundaries.

    Parameters
    ----------
    mask : np.ndarray
        Integer label image where 0 is background.
    seed_centroids : dict, optional
        Mapping of ``label -> (row, col)`` centroid to use as watershed
        seeds.  If ``None``, centroids are computed from ``mask`` itself
        via ``regionprops``.  Providing external centroids (e.g. from the
        nuclear mask) avoids using corrupted centroids from a flattened
        overlapping mask where region shapes are already distorted.

    Returns
    -------
    np.ndarray
        Re-partitioned label mask with the same dtype as *mask*.
    """
    props = label_props_dict(mask)
    if not props:
        return mask

    seeds = np.zeros_like(mask, dtype=mask.dtype)
    for lab in props:
        if seed_centroids and lab in seed_centroids:
            r, c = seed_centroids[lab]
        else:
            r, c = props[lab]["centroid"]
        r = max(0, min(int(round(r)), mask.shape[0] - 1))
        c = max(0, min(int(round(c)), mask.shape[1] - 1))
        seeds[r, c] = lab

    occupied = mask > 0
    dist = ndi.distance_transform_edt(seeds == 0)
    result = watershed(dist, markers=seeds, mask=occupied).astype(mask.dtype)
    return result


def main() -> int:
    """Entry point: orchestrate the full cell measurement pipeline.

    **Pipeline stages:**

    1. **Parse arguments** and validate constraints.
    2. **Load inputs** — nuclear mask, whole-cell mask, multi-channel image.
       Optionally downsample all three arrays by the same factor.
    3. **Match cells** — pair nuclear and whole-cell labels via centroid
       proximity (``match_cells``); synthesize boundaries for unmatched
       nuclei via watershed.  In ``--skip-nuclear-mask`` mode, whole-cell
       labels are used directly.
    4. **Measure** — for each cell, extract a bounding-box crop and compute
       shape metrics + intensity statistics from the **raw raster mask**
       (``feature_for_cell``, possibly parallelised across
       ``--threads`` workers).  Results are streamed to a temporary JSONL
       file to limit peak memory.
    5. **Free arrays** — release the large image and mask arrays to reclaim
       memory before post-processing.
    6. **Overlap clipping** — resolve overlapping cell polygons so no two
       cells share area (``constrain_cell_overlaps``).  Measurements are
       **not** recomputed; they reflect the raw segmentation masks.
    7. **Neighbourhood features** — optionally aggregate each cell's
       measurements across its k nearest neighbours.
    8. **Export** — write the GeoJSON FeatureCollection (with an
       annotation feature for the whole image extent) and optionally a
       rasterized label mask TIFF from the clipped polygons.

    Returns
    -------
    int
        Exit code (0 on success).
    """
    args = parse_args()
    # ===================================================================
    # Stage 1: Validate arguments
    # ===================================================================
    if args.tile_size <= 0:
        raise ValueError("tile-size must be > 0")
    if args.tile_overlap < 0:
        raise ValueError("tile-overlap must be >= 0")
    if args.pixel_size_microns <= 0:
        raise ValueError("pixel-size-microns must be > 0")
    if args.dist_threshold <= 0:
        raise ValueError("dist-threshold must be > 0")
    if args.downsample_factor <= 0:
        raise ValueError("downsample-factor must be > 0")
    step = int(round(args.downsample_factor))
    if args.downsample_factor > 1.0 and step < 2:
        print(
            f"Warning: downsample-factor {args.downsample_factor} rounds to step={step}; no downsampling will be applied"
        )
    if tile_flags_explicit(sys.argv[1:]):
        print("Warning: --tile-size/--tile-overlap are parsed for compatibility but not used in this Python implementation")

    # ===================================================================
    # Stage 2: Load inputs and optionally downsample
    # ===================================================================
    whole = load_label_mask(args.whole_cell_mask)
    if args.skip_nuclear_mask:
        nuc = None
    else:
        nuc = load_label_mask(args.nuclear_mask)
    img_cyx, ch_names = load_image(args.tiff_file)

    if nuc is not None:
        img_cyx, nuc, whole = maybe_downsample(img_cyx, nuc, whole, args.downsample_factor)
        if whole.shape != nuc.shape:
            raise ValueError(f"Mask shapes differ: whole={whole.shape}, nuclear={nuc.shape}")
    else:
        img_cyx, _, whole = maybe_downsample(img_cyx, whole, whole, args.downsample_factor)
    if img_cyx.shape[1:] != whole.shape:
        raise ValueError(f"Image shape {img_cyx.shape[1:]} does not match mask shape {whole.shape}")

    print(f"Loaded whole cell mask: {whole.shape}")
    if nuc is not None:
        print(f"Loaded nuclear mask: {nuc.shape}")
    else:
        print("Nuclear mask: skipped (--skip-nuclear-mask)")

    if args.skip_nuclear_mask:
        print("--skip-nuclear-mask: using whole-cell mask only (no compartmental measurements)")
        # Use whole-cell labels directly; broadcast stub avoids allocating a full-size zeros array (~6 GiB).
        cell_labels = whole
        nuc_labels = np.broadcast_to(np.int64(0), whole.shape)
        unique_cells = [int(x) for x in np.unique(cell_labels) if x > 0]
        # Build bbox map from the whole-cell label mask
        from skimage.measure import regionprops as _rp
        bbox_map = {}
        for r in _rp(cell_labels):
            bbox_map[r.label] = (slice(r.bbox[0], r.bbox[2]), slice(r.bbox[1], r.bbox[3]))
        records_by_id = {cid: CellRecord(cell_id=cid, cell_label=cid, nucleus_label=None) for cid in unique_cells}
        match_stats = {
            'nucleus_count': 0, 'whole_cell_count': len(unique_cells),
            'matched_cells': 0, 'unmatched_whole_cells': len(unique_cells),
            'dropped_synth_cells': 0,
        }
        image_shape = whole.shape
    else:
        cell_labels, nuc_labels, records, match_stats, bbox_map = match_cells(
            nuc, whole, args.dist_threshold, args.estimate_cell_boundary_dist
        )
        records_by_id = {r.cell_id: r for r in records}
        image_shape = whole.shape
        del nuc, whole
        gc.collect()

    unique_cells = [int(x) for x in np.unique(cell_labels) if x > 0]
    bbox_ids = set(bbox_map.keys())
    if len(unique_cells) != len(bbox_ids) or set(unique_cells) != bbox_ids:
        raise RuntimeError(
            "Internal mismatch between labeled cells and bbox map keys; "
            f"labels={len(unique_cells)}, bboxes={len(bbox_ids)}"
        )
    print(f"Total path objects: {len(unique_cells)}")
    print(
        "Matching summary: "
        f"nuclei={match_stats['nucleus_count']}, whole_cells={match_stats['whole_cell_count']}, "
        f"matched={match_stats['matched_cells']}, unmatched_whole={match_stats['unmatched_whole_cells']}, "
        f"dropped_estimated={match_stats['dropped_synth_cells']}"
    )

    percentiles = parse_csv_numbers(args.percentiles, cast=float)

    if percentiles:
        print(f"Will add intensity percentiles: {percentiles}")
    if args.erosion_steps:
        print(f"Will add erosion bin measurements (5 equal-area bins)")
    expand_20um_px = max(1, int(round(20.0 / args.pixel_size_microns)))
    if args.expansion_steps:
        print(f"Will add expansion bin measurements (20 µm = {expand_20um_px} px, 5 bins)")
    if args.environment_expansion:
        print(f"Will add environment measurements (20 µm = {expand_20um_px} px dilation)")

    # ===================================================================
    # Stage 4: Per-cell measurement (parallelisable)
    # Stream features to a temp JSONL file to limit peak memory.
    # ===================================================================
    out_path = Path(args.output_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    _tmp_path = str(out_path) + '.features.jsonl.tmp'
    _tmp_f = open(_tmp_path, 'w', encoding='utf-8')
    feature_count = 0
    total = len(unique_cells)

    # Expose large arrays via module globals so fork()-based worker processes
    # inherit them via copy-on-write instead of pickling per-task copies.
    # Workers only read from these globals, so no CoW page faults occur.
    global _GLOBAL_IMG, _GLOBAL_CELL, _GLOBAL_NUC, _GLOBAL_SKIP_NUC  # declared here; cleared again in Stage 5
    _GLOBAL_IMG = img_cyx
    _GLOBAL_CELL = cell_labels
    _GLOBAL_NUC = None if args.skip_nuclear_mask else nuc_labels
    _GLOBAL_SKIP_NUC = args.skip_nuclear_mask

    h_img, w_img = image_shape
    task_iter = iter_tasks_coords(
        unique_cells,
        bbox_map,
        records_by_id,
        h_img,
        w_img,
        args,
        ch_names,
        percentiles,
        args.erosion_steps,
        args.expansion_steps,
        args.environment_expansion,
    )

    if args.threads > 1:
        # --- Multi-process execution: bounded work queue to limit memory ---
        with ProcessPoolExecutor(max_workers=args.threads) as ex:
            max_inflight = max(1, args.threads * 4)
            future_map = {}
            for _ in range(max_inflight):
                try:
                    task = next(task_iter)
                except StopIteration:
                    break
                future_map[ex.submit(_feature_for_cell_global, *task)] = task[0]

            done = 0
            while future_map:
                completed, _ = wait(future_map.keys(), return_when=FIRST_COMPLETED)
                for fut in completed:
                    cid = future_map.pop(fut)
                    try:
                        feat = fut.result()
                    except Exception as exc:
                        print(f"Warning: cell {cid} failed: {exc}", file=sys.stderr)
                        feat = None

                    if feat is not None:
                        json.dump(feat, _tmp_f, separators=(',', ':'))
                        _tmp_f.write('\n')
                        feature_count += 1
                    done += 1
                    if done % 1000 == 0 or done == total:
                        print(f"Progress: {done}/{total} cells ({int(done * 100.0 / max(total, 1))}%)")

                while len(future_map) < max_inflight:
                    try:
                        task = next(task_iter)
                    except StopIteration:
                        break
                    future_map[ex.submit(_feature_for_cell_global, *task)] = task[0]
    else:
        # --- Single-threaded execution: simpler, easier to debug ---
        for i, t in enumerate(task_iter, 1):
            feat = _feature_for_cell_global(*t)
            if feat is not None:
                json.dump(feat, _tmp_f, separators=(',', ':'))
                _tmp_f.write('\n')
                feature_count += 1
            if i % 1000 == 0 or i == total:
                print(f"Progress: {i}/{total} cells ({int(i * 100.0 / max(total, 1))}%)")

    _tmp_f.close()

    # ===================================================================
    # Stage 5: Free large arrays to reclaim memory for post-processing
    # ===================================================================
    del task_iter, img_cyx, cell_labels, nuc_labels, bbox_map, records_by_id
    # Clear module globals so the arrays can be garbage collected.
    _GLOBAL_IMG = _GLOBAL_CELL = _GLOBAL_NUC = None
    gc.collect()
    print(f"Wrote {feature_count} features to temp file; freed image/mask arrays for post-processing")

    # ===================================================================
    # Stage 6: Load features and resolve geometry overlaps
    # Note: measurements are NOT recomputed after clipping — they reflect
    # the raw segmentation masks, which is intentional (see docstrings).
    # ===================================================================
    print("Loading features from temp file...")
    features = []
    with open(_tmp_path, 'r', encoding='utf-8') as _tmp_r:
        for line in _tmp_r:
            line = line.strip()
            if line:
                features.append(json.loads(line))
    os.unlink(_tmp_path)
    print(f"Loaded {len(features)} features")

    # Keep output deterministic regardless of parallel execution completion order.
    features.sort(key=lambda f: f["properties"].get("id", -1))

    # Resolve overlapping cell geometries (equivalent to QuPath CellTools.constrainCellOverlaps)
    features = constrain_cell_overlaps(features)

    # ===================================================================
    # Stage 7: Neighbourhood feature aggregation (optional)
    # Uses post-clipping centroids for spatial proximity but pre-clipping
    # measurement values for the aggregated statistics.
    # ===================================================================
    if args.neighbors and args.neighbors > 0:
        print(f"Computing neighborhood features (k={args.neighbors}, max 20 µm = {20.0 / args.pixel_size_microns:.1f} px)...")
        add_neighborhood_features(features, args.neighbors, pixel_size_microns=args.pixel_size_microns)
        print(f"Neighborhood features added for {len(features)} cells")

    # ===================================================================
    # Stage 8: Assemble and export GeoJSON FeatureCollection
    # ===================================================================
    # Top-level annotation feature for whole image extent — acts as a
    # bounding region when the GeoJSON is loaded into QuPath.
    h, w = image_shape
    annotation = {
        "type": "Feature",
        "id": "annotation-whole-image",
        "geometry": mapping(Polygon([(0, 0), (w, 0), (w, h), (0, h), (0, 0)])),
        "properties": {
            "objectType": "annotation",
            "type": "annotation",
            "name": "whole_image",
        },
    }

    out = {
        "type": "FeatureCollection",
        "features": [annotation] + features,
    }

    if args.gzip:
        if not str(out_path).endswith(".gz"):
            out_path = Path(str(out_path) + ".gz")
        with gzip.open(out_path, "wt", encoding="utf-8") as f:
            if args.pretty_json:
                json.dump(out, f, indent=2)
            else:
                json.dump(out, f, separators=(",", ":"))
    else:
        with out_path.open("w", encoding="utf-8") as f:
            if args.pretty_json:
                json.dump(out, f, indent=2)
            else:
                json.dump(out, f, separators=(",", ":"))

    print(f"Exported to GeoJSON: {out_path}")

    # Optionally rasterize the clipped polygons back to a label mask TIFF.
    # This mask reflects the POST-clipping geometry (non-overlapping), unlike
    # the measurements which reflect the pre-clipping raw segmentation masks.
    if args.output_mask:
        mask_out = rasterize_features_to_mask(features, h, w)
        mask_path = Path(args.output_mask)
        mask_path.parent.mkdir(parents=True, exist_ok=True)
        tifffile.imwrite(str(mask_path), mask_out)
        print(f"Exported label mask: {mask_path} ({mask_out.dtype}, {len(np.unique(mask_out)) - 1} labels)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
