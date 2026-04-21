#!/usr/bin/env python3
"""
Smooth label segmentation masks to reduce polygon complexity.

Two methods are supported, both parallelised via shared memory for scalability
on large masks (e.g. 40k x 40k with 100k+ labels):

  morphological (default)
    Morphological close + open per label using a disk structuring element.

  shapely
    Converts each label's boundary to a Shapely polygon and runs
    shapely.simplify (Douglas-Peucker) to reduce vertex count directly.
    This is the most direct way to reduce polygon complexity before GeoJSON
    export and can prevent StackOverflowError in QuPath.

Both methods write a new integer label TIFF of the same dtype as the input.
"""

import argparse
import sys
import numpy as np
import tifffile
from multiprocessing import shared_memory
from concurrent.futures import ProcessPoolExecutor
from scipy.ndimage import find_objects
from skimage.morphology import disk, binary_closing, binary_opening
from skimage.measure import find_contours
from skimage.draw import polygon as draw_polygon
from shapely.geometry import Polygon
from shapely.ops import unary_union


def enforce_min_area(mask, min_area):
    """Remove labels whose final pixel area is below min_area."""
    if min_area <= 0:
        return mask

    labels, counts = np.unique(mask, return_counts=True)
    small = labels[(labels != 0) & (counts < min_area)]
    if small.size == 0:
        return mask

    out = mask.copy()
    out[np.isin(out, small)] = 0
    return out


# ── Worker (must be module-level for pickling) ────────────────────────────────

def _process_label(args):
    """
    Process a single label: extract crop → morphological close+open → return pixels.
    Operates in crop-local space then offsets to global coordinates.
    """
    (label_id, slc_starts, slc_stops, shm_name,
     mask_shape, mask_dtype, kernel_size, min_area) = args

    shm = shared_memory.SharedMemory(name=shm_name)
    label_mask = np.ndarray(mask_shape, dtype=mask_dtype, buffer=shm.buf)

    try:
        # Pad bounding box by kernel_size so the disk structuring element isn't
        # clipped at the crop edge for labels near the image interior boundaries.
        H, W = mask_shape[0], mask_shape[1]
        pad_starts = tuple(max(0, s - kernel_size) for s in slc_starts)
        pad_stops  = tuple(min(d, e + kernel_size) for e, d in zip(slc_stops, (H, W)))
        slc = tuple(slice(s, e) for s, e in zip(pad_starts, pad_stops))
        crop = label_mask[slc] == label_id

        if not np.any(crop):
            return label_id, None, None

        selem = disk(kernel_size)

        # Close first (fills small holes, connects nearby regions)
        binary = binary_closing(crop, selem)
        # Then open (removes small protrusions, smooths boundary)
        binary = binary_opening(binary, selem)

        if not np.any(binary):
            return label_id, None, None

        # Get local pixel coordinates and offset to global using padded origin.
        local_rr, local_cc = np.where(binary)
        minr, minc = pad_starts
        rr = local_rr + minr
        cc = local_cc + minc

        # Clip to mask bounds — guards against kernel expanding beyond crop edge
        valid = (
            (rr >= 0) & (rr < mask_shape[0]) &
            (cc >= 0) & (cc < mask_shape[1])
        )

        # Enforce minimum area on final rasterized pixels.
        if np.count_nonzero(valid) < min_area:
            return label_id, None, None

        return label_id, rr[valid], cc[valid]

    except Exception as e:
        print(f"  [warn] label {label_id} failed: {e}", flush=True)
        # Fallback: return original pixels unchanged
        orig_rr, orig_cc = np.where(label_mask == label_id)
        if len(orig_rr) < min_area:
            return label_id, None, None
        return label_id, orig_rr, orig_cc

    finally:
        shm.close()  # Never unlink in workers — only the creator does that


# ── Main smoothing function ───────────────────────────────────────────────────

def smooth_label_morphological_parallel(
    label_mask, kernel_size=2, min_area=0, n_workers=8, chunksize=200
):
    """
    Parallelised morphological label smoothing using shared memory.

    Parameters
    ----------
    label_mask : np.ndarray
        Integer label image where 0 is background.
    kernel_size : int
        Radius of the disk structuring element (default: 2).
    min_area : float
        Minimum area in pixels^2 to retain after smoothing.
        Labels with fewer pixels than this threshold are removed.
    n_workers : int
        Number of parallel worker processes.
    chunksize : int
        Number of labels per executor dispatch chunk.
    """
    print(f"Setting up shared memory for mask {label_mask.shape} {label_mask.dtype}...")
    shm = shared_memory.SharedMemory(create=True, size=label_mask.nbytes)
    shared_arr = np.ndarray(label_mask.shape, dtype=label_mask.dtype, buffer=shm.buf)
    np.copyto(shared_arr, label_mask)

    print("Finding label bounding boxes...")
    slices = find_objects(label_mask)
    n_labels = sum(1 for s in slices if s is not None)
    print(f"Processing {n_labels} labels with {n_workers} workers, chunksize={chunksize}...")

    # Serialise slices as (starts, stops) tuples — slice objects aren't picklable
    work = []
    for label_id, slc in enumerate(slices, start=1):
        if slc is None:
            continue
        slc_starts = tuple(s.start for s in slc)
        slc_stops  = tuple(s.stop  for s in slc)
        work.append((
            label_id, slc_starts, slc_stops,
            shm.name, label_mask.shape, label_mask.dtype,
            kernel_size, min_area
        ))

    smoothed = np.zeros_like(label_mask)
    done = 0
    warned = 0

    try:
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            for label_id, rr, cc in executor.map(
                _process_label, work, chunksize=chunksize
            ):
                if rr is not None and len(rr) > 0:
                    # First-write wins for any overlap between adjacent labels
                    free = smoothed[rr, cc] == 0
                    smoothed[rr[free], cc[free]] = label_id
                else:
                    warned += 1
                done += 1
                if done % 5000 == 0:
                    print(f"  {done}/{n_labels} labels done...", flush=True)
    finally:
        shm.close()
        shm.unlink()  # Only unlink in the creator process

    if warned > 0:
        print(f"  [warn] {warned} labels produced no output pixels (removed by morphological op or too small)")

    smoothed = enforce_min_area(smoothed, min_area)

    return smoothed


# ── Shapely worker (must be module-level for pickling) ────────────────────────

def _process_label_shapely(args):
    """
    Process a single label using Shapely polygon simplification.

    Pipeline: binary crop → skimage contour → Shapely Polygon →
              shapely.simplify(tolerance) → skimage rasterize → global pixels.
    """
    (label_id, slc_starts, slc_stops, shm_name,
     mask_shape, mask_dtype, tolerance, min_area) = args

    shm = shared_memory.SharedMemory(name=shm_name)
    label_mask = np.ndarray(mask_shape, dtype=mask_dtype, buffer=shm.buf)

    try:
        slc = tuple(slice(s, e) for s, e in zip(slc_starts, slc_stops))
        crop = (label_mask[slc] == label_id).astype(np.uint8)

        if not np.any(crop):
            return label_id, None, None

        # Pad so contours touching the crop border are closed cleanly
        padded = np.pad(crop, 1, mode='constant', constant_values=0)
        contours = find_contours(padded, 0.5)

        if not contours:
            return label_id, None, None

        # Build a Shapely polygon for every contour (undo the 1-pixel pad offset).
        # Using unary_union handles disconnected fragments and fills holes, both of
        # which are acceptable for a smoothing operation on segmentation masks.
        polys = []
        for c in contours:
            c_shifted = c - 1.0
            if len(c_shifted) < 3:
                continue
            p = Polygon(zip(c_shifted[:, 1], c_shifted[:, 0]))
            if not p.is_valid:
                p = p.buffer(0)
            if not p.is_empty and p.area > 0:
                polys.append(p)

        if not polys:
            return label_id, None, None

        poly = unary_union(polys)
        poly = poly.simplify(tolerance, preserve_topology=True)

        if poly.is_empty:
            return label_id, None, None

        if poly.area < min_area:
            return label_id, None, None

        # Collect all parts — simplification can produce a MultiPolygon.
        if poly.geom_type == 'MultiPolygon':
            parts = list(poly.geoms)
        elif poly.geom_type == 'Polygon':
            parts = [poly]
        else:
            return label_id, None, None

        crop_h = slc_stops[0] - slc_starts[0]
        crop_w = slc_stops[1] - slc_starts[1]
        local_mask = np.zeros((crop_h, crop_w), dtype=bool)

        for part in parts:
            if part.is_empty:
                continue
            ext = np.array(part.exterior.coords)
            rr_e, cc_e = draw_polygon(ext[:, 1], ext[:, 0], shape=(crop_h, crop_w))
            local_mask[rr_e, cc_e] = True

        local_rr, local_cc = np.where(local_mask)

        # Offset to global coordinates
        rr = local_rr + slc_starts[0]
        cc = local_cc + slc_starts[1]

        # Clip to mask bounds
        valid = (
            (rr >= 0) & (rr < mask_shape[0]) &
            (cc >= 0) & (cc < mask_shape[1])
        )

        if not np.any(valid):
            return label_id, None, None

        # Enforce threshold on final rasterized pixel area so behavior matches
        # downstream measured mask areas.
        if np.count_nonzero(valid) < min_area:
            return label_id, None, None

        return label_id, rr[valid], cc[valid]

    except Exception as e:
        print(f"  [warn] label {label_id} failed: {e}", flush=True)
        # Fallback: return original pixels unchanged
        orig_rr, orig_cc = np.where(label_mask == label_id)
        if len(orig_rr) < min_area:
            return label_id, None, None
        return label_id, orig_rr, orig_cc

    finally:
        shm.close()


# ── Shapely smoothing driver ──────────────────────────────────────────────────

def smooth_label_shapely_parallel(
    label_mask, tolerance=1.0, min_area=0, n_workers=8, chunksize=200
):
    """
    Parallelised Shapely polygon simplification using shared memory.

    Each label boundary is extracted as a Shapely Polygon, simplified with the
    Douglas-Peucker algorithm (``shapely.simplify``), then rasterized back to
    pixel coordinates.  This directly reduces polygon vertex count, which is the
    root cause of StackOverflowError in QuPath's GeoJSON export.

    Parameters
    ----------
    label_mask : np.ndarray
        Integer label image where 0 is background.
    tolerance : float
        Douglas-Peucker simplification tolerance in pixels (default: 1.0).
        Higher values produce simpler polygons but less accurate boundaries.
    min_area : float
        Minimum polygon area in pixels² to retain after simplification.
        Labels whose simplified polygon falls below this are removed (default: 0).
    n_workers : int
        Number of parallel worker processes.
    chunksize : int
        Number of labels per executor dispatch chunk.
    """
    print(f"Setting up shared memory for mask {label_mask.shape} {label_mask.dtype}...")
    shm = shared_memory.SharedMemory(create=True, size=label_mask.nbytes)
    shared_arr = np.ndarray(label_mask.shape, dtype=label_mask.dtype, buffer=shm.buf)
    np.copyto(shared_arr, label_mask)

    print("Finding label bounding boxes...")
    slices = find_objects(label_mask)
    n_labels = sum(1 for s in slices if s is not None)
    print(f"Processing {n_labels} labels with {n_workers} workers, chunksize={chunksize}...")

    work = []
    for label_id, slc in enumerate(slices, start=1):
        if slc is None:
            continue
        slc_starts = tuple(s.start for s in slc)
        slc_stops  = tuple(s.stop  for s in slc)
        work.append((
            label_id, slc_starts, slc_stops,
            shm.name, label_mask.shape, label_mask.dtype,
            tolerance, min_area
        ))

    smoothed = np.zeros_like(label_mask)
    done = 0
    warned = 0

    try:
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            for label_id, rr, cc in executor.map(
                _process_label_shapely, work, chunksize=chunksize
            ):
                if rr is not None and len(rr) > 0:
                    free = smoothed[rr, cc] == 0
                    smoothed[rr[free], cc[free]] = label_id
                else:
                    warned += 1
                done += 1
                if done % 5000 == 0:
                    print(f"  {done}/{n_labels} labels done...", flush=True)
    finally:
        shm.close()
        shm.unlink()

    if warned > 0:
        print(f"  [warn] {warned} labels produced no output pixels (removed by simplification or too small)")

    smoothed = enforce_min_area(smoothed, min_area)

    return smoothed


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Smooth label segmentation masks to reduce polygon complexity."
    )
    parser.add_argument("input_mask",  help="Input label mask TIFF.")
    parser.add_argument("output_mask", help="Output smoothed mask TIFF.")

    # --- method selection ---
    parser.add_argument(
        "--method", choices=["morphological", "shapely"],
        default="morphological",
        help="Smoothing method: 'morphological' (close+open, default) or "
             "'shapely' (Douglas-Peucker polygon simplification)."
    )

    # --- morphological options ---
    parser.add_argument(
        "--kernel-size", type=int, default=2,
        help="[morphological] Radius of the disk structuring element for "
             "morphological close+open (default: 2)."
    )

    # --- shapely options ---
    parser.add_argument(
        "--tolerance", type=float, default=1.0,
        help="[shapely] Douglas-Peucker simplification tolerance in pixels "
             "(default: 1.0). Higher values produce simpler polygons."
    )
    parser.add_argument(
        "--min-area", type=float, default=0.0,
           help="Minimum area threshold in pixels\u00b2. Labels are dropped if the "
               "final rasterized pixel area falls below this value. In shapely "
               "mode, simplified polygon area is also checked (default: 0)."
    )

    # --- shared options ---
    parser.add_argument(
        "--n-workers", type=int, default=8,
        help="Number of parallel worker processes (default: 8). "
             "Match to --cpus-per-task in your SLURM script."
    )
    parser.add_argument(
        "--chunksize", type=int, default=200,
        help="Labels per worker chunk (default: 200)."
    )
    args = parser.parse_args()

    print(f"Reading mask: {args.input_mask}")
    mask = tifffile.imread(args.input_mask)
    original_dtype = mask.dtype
    print(f"Mask shape: {mask.shape}, dtype: {original_dtype}")

    n_labels = len(np.unique(mask)) - 1  # exclude background
    print(f"Number of labels: {n_labels}")

    if n_labels == 0:
        print("No labels found — writing input unchanged.")
        tifffile.imwrite(args.output_mask, mask)
        sys.exit(0)

    print(f"Using smoothing method: {args.method}")

    if args.method == "morphological":
        smoothed = smooth_label_morphological_parallel(
            mask,
            kernel_size=args.kernel_size,
            min_area=args.min_area,
            n_workers=args.n_workers,
            chunksize=args.chunksize,
        )
    else:  # shapely
        smoothed = smooth_label_shapely_parallel(
            mask,
            tolerance=args.tolerance,
            min_area=args.min_area,
            n_workers=args.n_workers,
            chunksize=args.chunksize,
        )

    # Preserve original dtype
    smoothed = smoothed.astype(original_dtype)

    n_after = len(np.unique(smoothed)) - 1
    print(f"Labels after smoothing: {n_after} (lost {n_labels - n_after})")

    print(f"Writing smoothed mask: {args.output_mask}")
    tifffile.imwrite(args.output_mask, smoothed)
    print("Done.")


if __name__ == "__main__":
    main()
