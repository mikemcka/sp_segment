#!/usr/bin/env python3
"""
KRONOS Embedding Extraction for Segmented Cells

This script extracts KRONOS foundation model embeddings for each segmented cell
in a multiplexed image. It takes the original TIFF image and whole-cell segmentation
mask, extracts cell-centered patches, matches marker channels to KRONOS marker
metadata (including fuzzy matching for unresolved markers), runs the KRONOS model,
and outputs per-cell embeddings as a CSV file.

Usage:
    python kronos_embeddings.py \
        --tiff image.tiff \
        --mask whole_cell_mask.tiff \
        --model-path /path/to/kronos_model.pt \
        --marker-metadata /path/to/marker_metadata.csv \
        --output embeddings.csv \
        [--patch-size 64] \
        [--batch-size 32] \
        [--max-value 65535] \
        [--marker-mapping '{"CD3e": "CD3E"}'] \
        [--geojson annotations.geojson] \
        [--merge-geojson]
"""

import argparse
import gc
import gzip
import json
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
import torch
import torch.utils.data
import tifffile

from difflib import SequenceMatcher
from scipy.ndimage import center_of_mass

try:
    from shapely.geometry import Polygon
    from rasterio import features
    GEOJSON_TO_MASK_AVAILABLE = True
except ImportError:
    GEOJSON_TO_MASK_AVAILABLE = False

warnings.filterwarnings("ignore")


def _mem_gb():
    """Return current process RSS memory in GB (Linux/macOS)."""
    try:
        import resource
        rss_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        # macOS returns bytes, Linux returns KB
        if sys.platform == "darwin":
            return rss_kb / (1024 ** 3)
        return rss_kb / (1024 ** 2)
    except Exception:
        return -1.0


def _log(msg):
    """Print a timestamped debug message with memory usage."""
    mem = _mem_gb()
    ts = time.strftime("%H:%M:%S")
    if mem >= 0:
        print(f"[{ts}] [RSS {mem:.1f} GB] {msg}", flush=True)
    else:
        print(f"[{ts}] {msg}", flush=True)


def get_channel_names(tif_path):
    """
    Extract channel/marker names from a TIFF file's metadata.

    Tries OME-XML metadata first, then ImageJ metadata, then falls back
    to numbered channel names.

    Args:
        tif_path: Path to the TIFF file.

    Returns:
        List of channel name strings.
    """
    _log(f"get_channel_names: opening {tif_path}")
    with tifffile.TiffFile(tif_path) as tif:
        # Get channel count from metadata without loading pixel data
        if tif.series:
            shape = tif.series[0].shape
            n_channels = shape[0] if len(shape) >= 3 else 1
        else:
            n_channels = len(tif.pages) if len(tif.pages) > 1 else 1
        _log(f"  n_channels={n_channels}, n_pages={len(tif.pages)}, is_ome={tif.is_ome}")

        channel_names = []

        # Try OME-XML metadata
        if tif.is_ome:
            try:
                from xml.etree import ElementTree as ET
                ome_xml = tif.ome_metadata
                root = ET.fromstring(ome_xml)
                ns = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
                channels = root.findall(".//ome:Channel", ns)
                if not channels:
                    channels = root.findall(".//Channel")
                channel_names = [ch.get("Name") or ch.get("ID") for ch in channels]
            except Exception:
                pass

        # Try ImageJ metadata
        if not channel_names and hasattr(tif, "imagej_metadata") and tif.imagej_metadata:
            if "Labels" in tif.imagej_metadata:
                channel_names = tif.imagej_metadata["Labels"]

        # Fallback to numbered channels
        if not channel_names:
            channel_names = [f"Channel_{i}" for i in range(n_channels)]

    return channel_names


def read_matched_channels(tiff_path, matched_channels):
    """
    Read only the specified channels from a TIFF file to avoid loading the
    entire image into memory. Falls back to full read + slice for single-page
    multi-channel TIFFs.

    Args:
        tiff_path: Path to the TIFF file.
        matched_channels: numpy array of channel indices to read.

    Returns:
        numpy array of shape (len(matched_channels), H, W).
    """
    t0 = time.time()
    with tifffile.TiffFile(tiff_path) as tif:
        n_pages = len(tif.pages)
        max_ch = int(max(matched_channels))
        _log(f"  TIFF has {n_pages} pages, need channels up to index {max_ch}")
        if n_pages > 1 and n_pages > max_ch:
            # Multi-page TIFF: read only needed pages (memory efficient)
            _log(f"  Using per-page reading strategy ({len(matched_channels)} pages)")
            first_page = tif.pages[int(matched_channels[0])].asarray()
            page_shape = first_page.shape
            page_dtype = first_page.dtype
            est_gb = (len(matched_channels) * np.prod(page_shape) * first_page.itemsize) / (1024**3)
            _log(f"  Page shape: {page_shape}, dtype: {page_dtype}, estimated total: {est_gb:.2f} GB")
            result = np.empty(
                (len(matched_channels), *page_shape),
                dtype=page_dtype,
            )
            result[0] = first_page
            del first_page
            for i, ch in enumerate(matched_channels[1:], 1):
                result[i] = tif.pages[int(ch)].asarray()
            _log(f"  Read {len(matched_channels)} pages in {time.time() - t0:.1f}s")
            return result
        else:
            # Single-page multi-channel or other layout: load full and slice
            _log(f"  Using full-load + slice strategy (n_pages={n_pages})")
            img = tif.asarray()
            if img.ndim == 2:
                img = img[np.newaxis, ...]
            _log(f"  Full image shape: {img.shape}, dtype: {img.dtype}, size: {img.nbytes / (1024**3):.2f} GB")
            result = img[matched_channels].copy()
            del img
            gc.collect()
            _log(f"  Sliced to {result.shape} in {time.time() - t0:.1f}s")
            return result


def match_markers(image_marker_names, kronos_metadata_df, user_mapping=None, top_suggestions=5):
    """
    Match image channel marker names to KRONOS marker metadata.

    Uses exact matching (case-insensitive), user-provided mappings, and
    fuzzy string matching (SequenceMatcher) for suggestions.

    Args:
        image_marker_names: List of marker names from the image channels.
        kronos_metadata_df: DataFrame of KRONOS marker metadata (marker_name, marker_id, marker_mean, marker_std).
        user_mapping: Optional dict mapping image marker names to KRONOS marker names.
        top_suggestions: Number of fuzzy-match suggestions to report for unmatched markers.

    Returns:
        Tuple of (matched_df, unmatched_report):
            - matched_df: DataFrame with columns [image_channel, marker_name, marker_id, marker_mean, marker_std]
              for all matched markers.
            - unmatched_report: String describing unmatched markers and suggestions.
    """
    # Build lookup by uppercase name
    kronos_metadata_df = kronos_metadata_df.copy()
    kronos_upper_lookup = {}
    for _, row in kronos_metadata_df.iterrows():
        kronos_upper_lookup[row["marker_name"].upper()] = row

    kronos_names_upper = list(kronos_upper_lookup.keys())

    matched_records = []
    unmatched_markers = []

    for ch_idx, marker_name in enumerate(image_marker_names):
        resolved_name = marker_name

        # Check user mapping first
        if user_mapping and marker_name in user_mapping:
            resolved_name = user_mapping[marker_name]

        upper_name = resolved_name.upper()

        if upper_name in kronos_upper_lookup:
            row = kronos_upper_lookup[upper_name]
            matched_records.append({
                "image_channel": ch_idx,
                "marker_name": marker_name,
                "kronos_marker_name": row["marker_name"],
                "marker_id": int(row["marker_id"]),
                "marker_mean": float(row["marker_mean"]),
                "marker_std": float(row["marker_std"]),
            })
        else:
            unmatched_markers.append((ch_idx, marker_name))

    # Build unmatched report with fuzzy suggestions
    report_lines = []
    if unmatched_markers:
        report_lines.append(f"WARNING: {len(unmatched_markers)} marker(s) could not be matched to KRONOS metadata:")
        report_lines.append(f"  (Matched {len(matched_records)} of {len(image_marker_names)} markers)")
        report_lines.append("")
        for ch_idx, marker_name in unmatched_markers:
            similarities = [
                (kn, SequenceMatcher(None, marker_name.upper(), kn).ratio())
                for kn in kronos_names_upper
            ]
            similarities.sort(key=lambda x: x[1], reverse=True)
            suggestions = [s[0] for s in similarities[:top_suggestions]]
            report_lines.append(f"  Channel {ch_idx}: '{marker_name}' -> suggestions: {suggestions}")
        report_lines.append("")
        report_lines.append("  To resolve, provide --marker-mapping as JSON, e.g.:")
        report_lines.append('    --marker-mapping \'{"MyMarker": "KRONOS_MARKER_NAME"}\'')

    matched_df = pd.DataFrame(matched_records)
    unmatched_report = "\n".join(report_lines)

    return matched_df, unmatched_report


def find_cell_centroids(mask):
    """
    Find cell IDs and centroids from a segmentation mask without extracting patches.

    Args:
        mask: numpy array of shape (H, W) - the whole-cell segmentation mask.

    Returns:
        Tuple of (cell_ids, centroids):
            - cell_ids: list of integer cell IDs
            - centroids: list of (y_center, x_center) tuples
    """
    t0 = time.time()
    _log(f"  Mask shape: {mask.shape}, dtype: {mask.dtype}, size: {mask.nbytes / (1024**3):.2f} GB")
    cell_ids_unique = np.unique(mask)
    cell_ids_unique = cell_ids_unique[cell_ids_unique > 0]  # skip background
    _log(f"  Found {len(cell_ids_unique)} unique cell labels (range: {cell_ids_unique.min()}-{cell_ids_unique.max()})" if len(cell_ids_unique) > 0 else "  No cells found in mask")

    if len(cell_ids_unique) == 0:
        return [], []

    # Use scipy for efficient vectorized centroid computation (single pass)
    centroids_array = center_of_mass(mask, labels=mask, index=cell_ids_unique)

    cell_ids = cell_ids_unique.tolist()
    centroids = [(int(round(c[0])), int(round(c[1]))) for c in centroids_array]
    _log(f"  Computed {len(cell_ids)} centroids in {time.time() - t0:.1f}s")

    return cell_ids, centroids


def load_kronos_model(model_path, config_path=None, device="cpu"):
    """
    Load the KRONOS model from a local checkpoint.

    Args:
        model_path: Path to the .pt model weights file.
        config_path: Path to config.json. If None, searches in the same directory as model_path.
        device: torch device string.

    Returns:
        Tuple of (model, precision, embedding_dim).
    """
    # Try to find config.json alongside the model
    if config_path is None:
        model_dir = os.path.dirname(model_path)
        candidate = os.path.join(model_dir, "config.json")
        if os.path.exists(candidate):
            config_path = candidate

    if config_path and os.path.exists(config_path):
        with open(config_path) as f:
            cfg = json.load(f)
    else:
        cfg = {"model_type": "vits16", "token_overlap": False}

    # Import kronos - it should be installed in the container/environment
    _log(f"  Importing kronos and loading checkpoint: {model_path}")
    _log(f"  Config path: {config_path}, config keys: {list(cfg.keys())}")
    from kronos import create_model_from_pretrained

    t0 = time.time()
    model, precision, embedding_dim = create_model_from_pretrained(
        checkpoint_path=model_path,
        cfg_path=config_path,
        cfg=cfg,
    )
    _log(f"  Model created in {time.time() - t0:.1f}s")

    model = model.to(device)
    model.eval()
    param_count = sum(p.numel() for p in model.parameters())
    _log(f"  Model on {device}: {param_count/1e6:.1f}M parameters, precision={precision}, embedding_dim={embedding_dim}")

    return model, precision, embedding_dim


def extract_embeddings(model, image, centroids, marker_ids, marker_means, marker_stds,
                       patch_size=64, max_value=65535.0, batch_size=32, num_workers=4,
                       device="cpu", precision=torch.float32):
    """
    Run KRONOS inference on cell patches to extract embeddings.

    Patches are extracted lazily from the image on-demand to avoid
    materializing all patches in memory at once.

    Args:
        model: Loaded KRONOS model.
        image: numpy array of shape (C, H, W) with the matched channels.
        centroids: list of (y_center, x_center) tuples for each cell.
        marker_ids: numpy array of marker IDs for each channel (length C).
        marker_means: numpy array of marker means for normalization (length C).
        marker_stds: numpy array of marker stds for normalization (length C).
        patch_size: Size of the square patch to extract around each cell centroid.
        max_value: Maximum intensity value for initial normalization.
        batch_size: Batch size for inference.
        num_workers: Number of DataLoader workers for parallel data loading.
        device: torch device string.
        precision: torch dtype for model input.

    Returns:
        Dict with keys:
            - 'patch_embeddings': numpy array of shape (N, embedding_dim)
            - 'marker_embeddings': numpy array of shape (N, num_markers, embedding_dim) or None
    """
    N = len(centroids)
    _log(f"extract_embeddings: N={N}, C={image.shape[0]}, patch_size={patch_size}, batch_size={batch_size}, num_workers={num_workers}")
    if N == 0:
        _log("  No centroids to process, returning empty")
        return {"patch_embeddings": np.array([]), "marker_embeddings": None}

    _log(f"  Image array: shape={image.shape}, dtype={image.dtype}, size={image.nbytes / (1024**3):.2f} GB")
    means = torch.tensor(marker_means, dtype=precision)
    stds = torch.tensor(marker_stds, dtype=precision)
    marker_ids_tensor = torch.tensor(marker_ids, dtype=torch.long)

    class _LazyPatchDataset(torch.utils.data.Dataset):
        """Extracts patches on-the-fly from the image to avoid storing all patches in RAM."""
        def __init__(self, image, centroids, patch_size, means, stds, max_value, precision):
            self.image = image  # (C, H, W) numpy array
            self.centroids = centroids
            self.patch_size = patch_size
            self.means = means
            self.stds = stds
            self.max_value = max_value
            self.precision = precision
            self.C, self.H, self.W = image.shape
            self.half = patch_size // 2

        def __len__(self):
            return len(self.centroids)

        def __getitem__(self, idx):
            y_center, x_center = self.centroids[idx]
            half = self.half

            y1 = y_center - half
            y2 = y_center + half
            x1 = x_center - half
            x2 = x_center + half

            patch = np.zeros((self.C, self.patch_size, self.patch_size), dtype=self.image.dtype)

            src_y1 = max(0, y1)
            src_y2 = min(self.H, y2)
            src_x1 = max(0, x1)
            src_x2 = min(self.W, x2)

            dst_y1 = src_y1 - y1
            dst_y2 = dst_y1 + (src_y2 - src_y1)
            dst_x1 = src_x1 - x1
            dst_x2 = dst_x1 + (src_x2 - src_x1)

            patch[:, dst_y1:dst_y2, dst_x1:dst_x2] = self.image[:, src_y1:src_y2, src_x1:src_x2]

            patch = torch.tensor(patch, dtype=self.precision)
            patch = patch / self.max_value
            patch = (patch - self.means[:, None, None]) / self.stds[:, None, None]
            return patch

    dataset = _LazyPatchDataset(image, centroids, patch_size, means, stds, max_value, precision)
    _log(f"  DataLoader: batch_size={batch_size}, num_workers={num_workers}, pin_memory={device != 'cpu'}")
    dataloader = torch.utils.data.DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=num_workers,
        pin_memory=(device != "cpu"),
        prefetch_factor=2 if num_workers > 0 else None,
    )

    all_patch_embeddings = []
    all_marker_embeddings = []
    processed = 0
    total_batches = len(dataloader)
    t_inference = time.time()

    for batch_idx, batch in enumerate(dataloader):
        batch = batch.to(device, non_blocking=True)
        batch_marker_ids = marker_ids_tensor.unsqueeze(0).expand(batch.shape[0], -1).to(device)

        with torch.no_grad():
            output = model(batch, marker_ids=batch_marker_ids)

        # Handle different output formats
        if isinstance(output, tuple):
            patch_emb = output[0].cpu().numpy()  # patch/cls embeddings
            marker_emb = output[1].cpu().numpy() if len(output) > 1 else None
        elif isinstance(output, dict):
            patch_emb = output.get("x_norm_clstoken", output.get("cls_token")).cpu().numpy()
            marker_emb = output.get("x_norm_patchtokens", output.get("patch_tokens"))
            if marker_emb is not None:
                marker_emb = marker_emb.cpu().numpy()
        else:
            patch_emb = output.cpu().numpy()
            marker_emb = None

        all_patch_embeddings.append(patch_emb)
        if marker_emb is not None:
            all_marker_embeddings.append(marker_emb)

        processed += batch.shape[0]
        if batch_idx == 0 or (batch_idx + 1) % max(1, total_batches // 10) == 0 or batch_idx == total_batches - 1:
            elapsed = time.time() - t_inference
            cells_per_sec = processed / elapsed if elapsed > 0 else 0
            _log(f"  Batch {batch_idx+1}/{total_batches}: {processed}/{N} cells ({cells_per_sec:.0f} cells/s)")

    total_inference_time = time.time() - t_inference
    _log(f"  Inference complete: {N} cells in {total_inference_time:.1f}s ({N/total_inference_time:.0f} cells/s)")

    result = {
        "patch_embeddings": np.concatenate(all_patch_embeddings, axis=0),
        "marker_embeddings": np.concatenate(all_marker_embeddings, axis=0) if all_marker_embeddings else None,
    }
    _log(f"  Embedding result: patch={result['patch_embeddings'].shape}, marker={'None' if result['marker_embeddings'] is None else result['marker_embeddings'].shape}")

    return result


def save_embeddings(cell_ids, centroids, embeddings, output_path, sample_id=None):
    """
    Save per-cell embeddings to a CSV file.

    Columns: cell_id, y_center, x_center, emb_0, emb_1, ..., emb_N

    Args:
        cell_ids: List of cell IDs.
        centroids: List of (y, x) centroid tuples.
        embeddings: numpy array of shape (num_cells, embedding_dim).
        output_path: Path for the output CSV file.
        sample_id: Optional sample identifier to include.
    """
    if len(cell_ids) == 0:
        # Write empty file with headers only
        with open(output_path, "w") as f:
            f.write("cell_id,y_center,x_center\n")
        _log(f"No cells found. Empty file written to {output_path}")
        return

    embedding_dim = embeddings.shape[1]
    emb_columns = [f"emb_{i}" for i in range(embedding_dim)]

    data = {
        "cell_id": cell_ids,
        "y_center": [c[0] for c in centroids],
        "x_center": [c[1] for c in centroids],
    }

    if sample_id:
        data["sample_id"] = [sample_id] * len(cell_ids)

    for i, col in enumerate(emb_columns):
        data[col] = embeddings[:, i]

    df = pd.DataFrame(data)
    df.to_csv(output_path, index=False)
    _log(f"Saved {len(cell_ids)} cell embeddings ({embedding_dim}D) to {output_path}")


def create_mask_from_geojson(geojson_path, width, height):
    """
    Create a segmentation mask from GeoJSON cell annotations.

    This ensures perfect 1:1 correspondence between mask labels and GeoJSON cells.
    Each cell in the GeoJSON is assigned a sequential label (1, 2, 3, ...).

    Args:
        geojson_path: Path to GeoJSON file with cell annotations
        width: Image width in pixels
        height: Image height in pixels

    Returns:
        Tuple of (mask, uuid_to_label_dict) where:
            - mask: numpy array of shape (height, width) with cell labels
            - uuid_to_label_dict: mapping from GeoJSON UUID to mask label
    """
    if not GEOJSON_TO_MASK_AVAILABLE:
        raise ImportError("shapely and rasterio are required for GeoJSON to mask conversion. "
                         "Install with: pip install shapely rasterio")

    _log(f"Creating mask from GeoJSON: {geojson_path}")

    open_func = gzip.open if str(geojson_path).endswith('.gz') else open
    with open_func(geojson_path, 'rt') as f:
        geojson = json.load(f)

    # Extract cell features and create label mapping
    cell_features = []
    uuid_to_label = {}
    label = 1

    for feature in geojson.get('features', []):
        if feature.get('properties', {}).get('objectType') == 'cell':
            cell_features.append(feature)
            uuid_to_label[feature['id']] = label
            label += 1

    _log(f"  Found {len(cell_features)} cells in GeoJSON")

    # Parse polygon coordinates
    def parse_polygon_coords(coords):
        # Flatten nested coordinate structures
        while len(coords) > 0 and isinstance(coords[0], list) and len(coords[0]) > 0 and isinstance(coords[0][0], list):
            coords = coords[0]
        return Polygon(coords)

    # Create shapes for rasterization
    shapes = []
    for feature in cell_features:
        coords = feature['geometry']['coordinates']
        poly = parse_polygon_coords(coords)
        label_val = uuid_to_label[feature['id']]
        shapes.append((poly, label_val))

    # Rasterize polygons into mask
    _log(f"  Rasterizing {len(shapes)} cells into {width}x{height} mask...")
    mask = features.rasterize(
        shapes,
        out_shape=(height, width),
        fill=0,
        dtype=np.uint32,
        all_touched=False
    )

    unique_labels = len(np.unique(mask)) - 1  # exclude background
    _log(f"  Created mask with {unique_labels} unique cell labels")

    return mask, uuid_to_label, geojson


def merge_embeddings_into_geojson(geojson_path, mask_path, cell_ids, centroids, embeddings, output_path,
                                  distance_threshold=5.0, use_geojson_mask=False, uuid_to_label=None,
                                  geojson_data=None):
    """
    Merge KRONOS embeddings into a GeoJSON file by matching via mask labels.

    Builds a mapping from GeoJSON cell UUID to mask label by extracting the mask
    value at each cell's centroid, then uses this to match embeddings to GeoJSON features.

    If use_geojson_mask=True, creates a new mask directly from the GeoJSON to ensure
    perfect 1:1 correspondence between mask labels and GeoJSON cells.

    Args:
        geojson_path: Path to the input GeoJSON from cellmeasurement.
        mask_path: Path to the whole-cell mask TIFF (or None if use_geojson_mask=True).
        cell_ids: List of cell IDs (mask labels) from embedding extraction.
        centroids: List of (y, x) centroid tuples from embedding extraction.
        embeddings: numpy array of shape (num_cells, embedding_dim).
        output_path: Path for the output merged GeoJSON.
        distance_threshold: Maximum distance (pixels) for centroid matching fallback.
        use_geojson_mask: If True, create mask from GeoJSON instead of loading from file.
        uuid_to_label: Dictionary mapping cell UUIDs to mask labels (used when use_geojson_mask=True).
        geojson_data: Pre-parsed GeoJSON dict to avoid reloading from disk.

    Returns:
        Tuple of (num_matched, num_unmatched, num_geojson_cells, method_counts).
    """
    if geojson_data is not None:
        geojson = geojson_data
        _log(f"  Using pre-parsed GeoJSON data (avoiding reload of {geojson_path})")
    else:
        open_func = gzip.open if str(geojson_path).endswith('.gz') else open
        with open_func(geojson_path, 'rt') as f:
            geojson = json.load(f)

    if use_geojson_mask:
        if uuid_to_label is None:
            # Create mask from GeoJSON when uuid_to_label was not pre-computed
            # First, get image dimensions from the original mask or from annotation
            if mask_path and os.path.exists(mask_path):
                orig_mask = tifffile.imread(mask_path)
                if orig_mask.ndim == 3:
                    orig_mask = orig_mask[0]
                mask_height, mask_width = orig_mask.shape
                del orig_mask
            else:
                # Get dimensions from annotation feature
                for feature in geojson.get('features', []):
                    if feature.get('properties', {}).get('objectType') == 'annotation':
                        coords = feature['geometry']['coordinates']
                        while len(coords) > 0 and isinstance(coords[0], list) and isinstance(coords[0][0], list):
                            coords = coords[0]
                        xs = [pt[0] for pt in coords]
                        ys = [pt[1] for pt in coords]
                        mask_width = int(max(xs))
                        mask_height = int(max(ys))
                        break

            # Create mask from GeoJSON
            _mask, uuid_to_label = create_mask_from_geojson(geojson_path, mask_width, mask_height)
            del _mask  # only need the uuid_to_label mapping

        # Build mask_label to embedding index mapping
        mask_label_to_emb_idx = {label: idx for idx, label in enumerate(cell_ids)}

        embedding_dim = embeddings.shape[1]
        emb_columns = [f"kronos_emb_{i}" for i in range(embedding_dim)]

        num_matched = 0
        num_unmatched = 0

        cell_features = [f for f in geojson.get('features', [])
                        if f.get('properties', {}).get('objectType') == 'cell']

        # Remove any existing KRONOS embeddings from a previous run
        cleared = 0
        for feature in cell_features:
            measurements = feature.get("properties", {}).get("measurements", {})
            old_keys = [k for k in measurements if k.startswith("kronos_")]
            if old_keys:
                for k in old_keys:
                    del measurements[k]
                cleared += 1
        if cleared > 0:
            _log(f"  Cleared existing KRONOS embeddings from {cleared} cells")

        for feature in cell_features:
            feature_id = feature["id"]

            if feature_id in uuid_to_label:
                mask_label = uuid_to_label[feature_id]
                if mask_label in mask_label_to_emb_idx:
                    emb_idx = mask_label_to_emb_idx[mask_label]
                    measurements = feature["properties"].setdefault("measurements", {})
                    measurements["kronos_cell_id"] = int(cell_ids[emb_idx])
                    for i, col_name in enumerate(emb_columns):
                        measurements[col_name] = float(embeddings[emb_idx, i])
                    num_matched += 1
                else:
                    num_unmatched += 1
            else:
                num_unmatched += 1

        _open = gzip.open if str(output_path).endswith('.gz') else open
        with _open(output_path, 'wt') as f:
            json.dump(geojson, f, separators=(',', ':'))

        method_counts = {"by_label": num_matched, "by_distance": 0}
        _log(f"  GeoJSON-derived mask matching: {num_matched}/{len(cell_features)} cells matched")

        if num_unmatched > 0:
            _log(f"  Warning: {num_unmatched} cells could not be matched")

        return num_matched, num_unmatched, len(cell_features), method_counts

    else:
        # Original mask-based matching logic
        # Load the whole-cell mask to extract labels at centroids
        mask = tifffile.imread(mask_path)
        if mask.ndim == 3:
            mask = mask[0]  # take first channel if multi-channel
        mask_height, mask_width = mask.shape

    # Extract centroids from GeoJSON cell features and map UUIDs to mask labels
    cell_features = []
    uuid_to_mask_label = {}
    unmapped_features = []

    for feature in geojson.get("features", []):
        props = feature.get("properties", {})
        if props.get("objectType") != "cell":
            continue
        cell_features.append(feature)

        # Compute centroid using Shoelace formula
        coords = feature["geometry"]["coordinates"]

        # Flatten nested coordinate structures to get to [x, y] pairs
        while len(coords) > 0 and isinstance(coords[0], list) and len(coords[0]) > 0 and isinstance(coords[0][0], list):
            coords = coords[0]  # unwrap one level

        # Now coords should be [[x, y], [x, y], ...]
        n = len(coords) - 1  # last point == first point (closed ring)
        A = 0.0
        cx = 0.0
        cy = 0.0

        for i in range(n):
            x0, y0 = coords[i][0], coords[i][1]
            x1, y1 = coords[i + 1][0], coords[i + 1][1]
            cross = x0 * y1 - x1 * y0
            A += cross
            cx += (x0 + x1) * cross
            cy += (y0 + y1) * cross
        A = A / 2.0
        if abs(A) > 1e-10:
            cx = cx / (6.0 * A)
            cy = cy / (6.0 * A)
        else:
            # Degenerate polygon, fall back to vertex average
            cx = sum(pt[0] for pt in coords) / len(coords)
            cy = sum(pt[1] for pt in coords) / len(coords)

        # Extract mask label at centroid (x, y) -> (col, row)
        cx_int = int(round(cx))
        cy_int = int(round(cy))

        if 0 <= cy_int < mask_height and 0 <= cx_int < mask_width:
            mask_label = int(mask[cy_int, cx_int])
            if mask_label > 0:  # skip background
                uuid_to_mask_label[feature["id"]] = mask_label
            else:
                unmapped_features.append((feature["id"], (cy_int, cx_int)))
        else:
            unmapped_features.append((feature["id"], (cy_int, cx_int)))

    if not cell_features or len(cell_ids) == 0:
        # Nothing to merge, write unchanged
        _open = gzip.open if str(output_path).endswith('.gz') else open
        with _open(output_path, 'wt') as f:
            json.dump(geojson, f, separators=(',', ':'))
        return 0, 0, len(cell_features), {}

    _log(f"  Mapped {len(uuid_to_mask_label)}/{len(cell_features)} GeoJSON cells to mask labels")
    if unmapped_features:
        _log(f"  Warning: {len(unmapped_features)} cells could not be mapped (centroid outside mask or on background)")

    # Build mask_label to embedding index mapping
    mask_label_to_emb_idx = {label: idx for idx, label in enumerate(cell_ids)}

    # Build mask_label to embedding index mapping
    mask_label_to_emb_idx = {label: idx for idx, label in enumerate(cell_ids)}

    embedding_dim = embeddings.shape[1]
    emb_columns = [f"kronos_emb_{i}" for i in range(embedding_dim)]

    num_matched = 0
    num_matched_by_label = 0
    num_matched_by_distance = 0
    num_unmatched = 0
    unmatched_distances = []

    # Remove any existing KRONOS embeddings from a previous run
    cleared = 0
    for feature in cell_features:
        measurements = feature.get("properties", {}).get("measurements", {})
        old_keys = [k for k in measurements if k.startswith("kronos_")]
        if old_keys:
            for k in old_keys:
                del measurements[k]
            cleared += 1
    if cleared > 0:
        _log(f"  Cleared existing KRONOS embeddings from {cleared} cells")

    # For each GeoJSON cell, try to match by mask label first
    for feature in cell_features:
        feature_id = feature["id"]
        matched = False

        # Try exact match via mask label
        if feature_id in uuid_to_mask_label:
            mask_label = uuid_to_mask_label[feature_id]
            if mask_label in mask_label_to_emb_idx:
                emb_idx = mask_label_to_emb_idx[mask_label]
                measurements = feature["properties"].setdefault("measurements", {})
                measurements["kronos_cell_id"] = int(cell_ids[emb_idx])
                for i, col_name in enumerate(emb_columns):
                    measurements[col_name] = float(embeddings[emb_idx, i])
                num_matched += 1
                num_matched_by_label += 1
                matched = True

        if not matched:
            num_unmatched += 1

    _open = gzip.open if str(output_path).endswith('.gz') else open
    with _open(output_path, 'wt') as f:
        json.dump(geojson, f, separators=(',', ':'))

    # Report matching statistics
    method_counts = {
        "by_label": num_matched_by_label,
        "by_distance": num_matched_by_distance
    }

    _log(f"  Matching statistics: {num_matched_by_label} by mask label, {num_matched_by_distance} by distance")

    # Report distance statistics for unmatched cells
    if unmatched_distances:
        unmatched_distances = np.array(unmatched_distances)
        _log(f"  Unmatched cell distance stats: min={unmatched_distances.min():.2f}, "
              f"median={np.median(unmatched_distances):.2f}, max={unmatched_distances.max():.2f} pixels")
        _log(f"  Consider increasing --distance-threshold (current: {distance_threshold}) if needed")

    return num_matched, num_unmatched, len(cell_features), method_counts


def main():
    parser = argparse.ArgumentParser(
        description="Extract KRONOS embeddings for segmented cells"
    )
    parser.add_argument("--tiff", required=True, help="Path to multiplexed TIFF image")
    parser.add_argument("--mask", required=True, help="Path to whole-cell segmentation mask TIFF")
    parser.add_argument("--model-path", required=True, help="Path to KRONOS model weights (.pt)")
    parser.add_argument("--config-path", default=None, help="Path to KRONOS config.json (auto-detected if not specified)")
    parser.add_argument("--marker-metadata", required=True, help="Path to KRONOS marker_metadata.csv")
    parser.add_argument("--output", required=True, help="Output CSV path for cell embeddings")
    parser.add_argument("--patch-size", type=int, default=64, help="Patch size for cell-centered crops (default: 64)")
    parser.add_argument("--batch-size", type=int, default=32, help="Batch size for model inference (default: 32)")
    parser.add_argument("--num-workers", type=int, default=4, help="Number of DataLoader workers for parallel data loading (default: 4)")
    parser.add_argument("--max-value", type=float, default=65535.0, help="Maximum intensity value for normalization (default: 65535)")
    parser.add_argument("--marker-mapping", default=None, help="JSON string or file path mapping image markers to KRONOS markers")
    parser.add_argument("--sample-id", default=None, help="Sample identifier to include in output")
    parser.add_argument("--geojson", default=None, help="Path to GeoJSON from cellmeasurement (required for --merge-geojson)")
    parser.add_argument("--merge-geojson", action="store_true", default=False,
                        help="Merge embeddings into GeoJSON as per-cell properties. Automatically creates mask from GeoJSON for perfect matching.")
    parser.add_argument("--output-geojson", default=None,
                        help="Explicit output path for merged GeoJSON. If not specified, derives from --output with _kronos suffix.")
    parser.add_argument("--distance-threshold", type=float, default=5.0,
                        help="Max pixel distance for centroid matching fallback (default: 5.0)")
    args = parser.parse_args()

    if args.merge_geojson and not args.geojson:
        parser.error("--merge-geojson requires --geojson to be specified")

    if args.merge_geojson and not GEOJSON_TO_MASK_AVAILABLE:
        parser.error("--merge-geojson requires shapely and rasterio for GeoJSON mask creation. Install with: pip install shapely rasterio")

    # Print all arguments for reproducibility
    _log("=" * 60)
    _log("KRONOS Embedding Extraction - Starting")
    _log("=" * 60)
    _log(f"  Arguments:")
    for k, v in vars(args).items():
        _log(f"    {k}: {v}")

    # Determine device
    device = "cuda" if torch.cuda.is_available() else "cpu"
    _log(f"Using device: {device}")
    if device == "cuda":
        _log(f"  CUDA device: {torch.cuda.get_device_name(0)}")
        _log(f"  CUDA memory: {torch.cuda.get_device_properties(0).total_memory / (1024**3):.1f} GB total")
    else:
        _log("  WARNING: No GPU detected, running on CPU")

    # Disable multiprocessing workers on CPU to avoid cleanup issues
    if device == "cpu" and args.num_workers > 0:
        _log(f"  Warning: Disabling num_workers (was {args.num_workers}) on CPU to avoid NFS/multiprocessing issues")
        args.num_workers = 0

    # Log file sizes for debugging
    for label, path in [("TIFF", args.tiff), ("Mask", args.mask), ("Marker metadata", args.marker_metadata),
                         ("GeoJSON", args.geojson), ("Marker mapping", args.marker_mapping)]:
        if path and os.path.isfile(path):
            size_mb = os.path.getsize(path) / (1024**2)
            _log(f"  File: {label} = {path} ({size_mb:.1f} MB)")

    # Load marker metadata
    _log(f"Loading KRONOS marker metadata from {args.marker_metadata}")
    kronos_metadata = pd.read_csv(args.marker_metadata)
    _log(f"  Metadata shape: {kronos_metadata.shape}, columns: {list(kronos_metadata.columns)}")

    # Parse user marker mapping
    user_mapping = None
    if args.marker_mapping:
        if os.path.isfile(args.marker_mapping):
            with open(args.marker_mapping) as f:
                user_mapping = json.load(f)
        else:
            try:
                user_mapping = json.loads(args.marker_mapping)
            except json.JSONDecodeError:
                _log(f"WARNING: Could not parse --marker-mapping as JSON: {args.marker_mapping}")

    # Get channel names from image metadata (reads metadata only, not pixel data)
    channel_names = get_channel_names(args.tiff)
    _log(f"  Found {len(channel_names)} channel names: {channel_names[:10]}{'...' if len(channel_names) > 10 else ''}")

    # Match markers
    matched_df, unmatched_report = match_markers(channel_names, kronos_metadata, user_mapping)
    if unmatched_report:
        _log(unmatched_report)
    _log(f"  Matched {len(matched_df)} markers to KRONOS metadata")

    if len(matched_df) == 0:
        _log("ERROR: No markers could be matched to KRONOS metadata. Exiting.")
        sys.exit(1)

    matched_channels = matched_df["image_channel"].values
    marker_ids = matched_df["marker_id"].values
    marker_means = matched_df["marker_mean"].values
    marker_stds = matched_df["marker_std"].values

    # Read only matched channels to avoid loading entire image into memory
    _log(f"Reading {len(matched_channels)} matched channels from: {args.tiff}")
    _log(f"  Channel indices: {matched_channels.tolist()}")
    image_matched = read_matched_channels(args.tiff, matched_channels)
    img_height, img_width = image_matched.shape[1], image_matched.shape[2]
    _log(f"  Matched image shape: {image_matched.shape} (H={img_height}, W={img_width}), size: {image_matched.nbytes / (1024**3):.2f} GB")

    # Read or create segmentation mask
    geojson_data = None  # keep parsed GeoJSON for merge step to avoid reloading
    if args.merge_geojson:
        # Create mask from GeoJSON for perfect matching when merging
        _log(f"Creating segmentation mask from GeoJSON: {args.geojson}")
        mask, uuid_to_label, geojson_data = create_mask_from_geojson(args.geojson, img_width, img_height)
        _log(f"  Mask shape: {mask.shape}, unique cells: {len(uuid_to_label)}, size: {mask.nbytes / (1024**3):.2f} GB")

        # Save the GeoJSON-derived mask for inspection/reuse
        prefix = args.sample_id if args.sample_id else "sample"
        geojson_mask_path = os.path.join(os.path.dirname(args.output), f"{prefix}_geojson_mask.tif")
        tifffile.imwrite(geojson_mask_path, mask.astype(np.uint32))
        _log(f"  Saved GeoJSON-derived mask to {geojson_mask_path}")
    else:
        # Read segmentation mask from file
        _log(f"Reading segmentation mask: {args.mask}")
        mask = tifffile.imread(args.mask)
        if mask.ndim == 3:
            mask = mask[0]  # take first channel if multi-channel
        _log(f"  Mask shape: {mask.shape}, dtype: {mask.dtype}, unique cells: {len(np.unique(mask)) - 1}, size: {mask.nbytes / (1024**3):.2f} GB")
        uuid_to_label = None

    # Extract cell centroids from mask (no patches stored in memory)
    _log(f"Finding cell centroids (patch_size={args.patch_size})...")
    cell_ids, centroids = find_cell_centroids(mask)
    _log(f"  Found {len(cell_ids)} cells")

    # Free mask memory - no longer needed (geojson merge recreates it if needed)
    _log("  Freeing mask memory")
    del mask
    gc.collect()

    if len(cell_ids) == 0:
        _log("WARNING: No cells found in segmentation mask.")
        save_embeddings([], [], np.array([]), args.output, args.sample_id)
        return

    # Load KRONOS model
    _log(f"Loading KRONOS model from {args.model_path}")
    model, precision, embedding_dim = load_kronos_model(args.model_path, args.config_path, device)
    _log(f"  Model loaded: precision={precision}, embedding_dim={embedding_dim}")
    if device == "cuda":
        _log(f"  GPU memory after model load: {torch.cuda.memory_allocated() / (1024**3):.2f} GB allocated, {torch.cuda.memory_reserved() / (1024**3):.2f} GB reserved")

    # Extract embeddings (patches extracted lazily to save memory)
    _log(f"Extracting embeddings (batch_size={args.batch_size}, num_workers={args.num_workers})...")
    results = extract_embeddings(
        model, image_matched, centroids, marker_ids, marker_means, marker_stds,
        patch_size=args.patch_size, max_value=args.max_value, batch_size=args.batch_size,
        num_workers=args.num_workers, device=device, precision=precision,
    )

    patch_embeddings = results["patch_embeddings"]
    _log(f"  Patch embeddings shape: {patch_embeddings.shape}")

    # Free model and image memory before GeoJSON merge
    _log("  Freeing model and image memory")
    del model, image_matched
    if device == "cuda":
        torch.cuda.empty_cache()
    gc.collect()
    _log(f"  Memory after cleanup")

    # Save embeddings
    _log(f"Saving embeddings to {args.output}")
    save_embeddings(cell_ids, centroids, patch_embeddings, args.output, args.sample_id)

    # Save marker matching report
    report_path = args.output.replace(".csv", "_marker_report.txt")
    with open(report_path, "w") as f:
        f.write("KRONOS Marker Matching Report\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Image: {args.tiff}\n")
        f.write(f"Total image channels: {len(channel_names)}\n")
        f.write(f"Matched markers: {len(matched_df)}\n\n")
        f.write("Matched markers:\n")
        for _, row in matched_df.iterrows():
            f.write(f"  Ch {row['image_channel']}: {row['marker_name']} -> {row['kronos_marker_name']} (ID={row['marker_id']})\n")
        if unmatched_report:
            f.write(f"\n{unmatched_report}\n")
    _log(f"  Marker matching report saved to {report_path}")

    # Optionally merge embeddings into GeoJSON
    if args.merge_geojson:
        if args.output_geojson:
            merged_geojson_path = args.output_geojson
        else:
            merged_geojson_path = args.output.replace("_kronos_embeddings.csv", "_kronos.geojson.gz")
            if not merged_geojson_path.endswith(".geojson.gz"):
                merged_geojson_path = args.output.rsplit(".", 1)[0] + "_kronos.geojson.gz"
        _log(f"Merging embeddings into GeoJSON: {args.geojson}")

        # Use GeoJSON-derived mask for perfect matching (pass pre-parsed data to avoid reloading 124GB+ file)
        num_matched, num_unmatched, num_geo_cells, method_counts = merge_embeddings_into_geojson(
            args.geojson, None, cell_ids, centroids, patch_embeddings,
            merged_geojson_path, distance_threshold=args.distance_threshold,
            use_geojson_mask=True, uuid_to_label=uuid_to_label,
            geojson_data=geojson_data
        )
        _log(f"  GeoJSON merge: {num_matched}/{num_geo_cells} cells matched, {num_unmatched} unmatched")
        _log(f"  Merged GeoJSON saved to {merged_geojson_path}")
        del geojson_data
        gc.collect()

    _log("=" * 60)
    _log("KRONOS Embedding Extraction - Done!")
    _log("=" * 60)


if __name__ == "__main__":
    main()
