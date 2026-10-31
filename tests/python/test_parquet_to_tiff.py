"""Unit tests for ``bin/parquet_to_tiff.py``: label dtype and pixel convention.

Label IDs are stored as pixel values, so a mask dtype too narrow for the largest
ID silently truncates the segmentation. This is worth pinning because there is
no error to notice: ``rasterio.features.rasterize`` does not raise on overflow,
it caps at the dtype maximum and every cell past that point vanishes. It only
bites on whole-slide images, which is exactly where nobody is watching.

Pure python -- no pipeline run required.
"""

import numpy as np
import rasterio.features
from shapely.geometry import box

import tifffile

from parquet_to_tiff import (
    LABEL_DTYPE,
    MASK_COMPRESSION,
    MASK_COMPRESSION_ARGS,
    rasterize_labels,
)

UINT16_MAX = np.iinfo(np.uint16).max  # 65535
UINT32_MAX = np.iinfo(np.uint32).max


def _label_shapes(n):
    """``n`` disjoint 1px cells with 1-based IDs, laid out in a 100-wide grid."""
    return ((box(i % 100, i // 100, i % 100 + 0.9, i // 100 + 0.9), i) for i in range(1, n + 1))


def _rasterize(n, dtype):
    height = (n // 100) + 2
    return rasterio.features.rasterize(
        _label_shapes(n), out_shape=(height, 110), fill=0, dtype=dtype, all_touched=True
    )


# ----------------------------------------------------------------------------
# the label dtype
# ----------------------------------------------------------------------------


def test_label_dtype_indexes_far_more_cells_than_any_real_slide():
    """uint32 is the wide fixed choice: 4.29e9 IDs, well past any real mask."""
    assert np.dtype(LABEL_DTYPE).name == "uint32"
    assert np.iinfo(LABEL_DTYPE).max == UINT32_MAX


# ----------------------------------------------------------------------------
# the truncation this guards against
# ----------------------------------------------------------------------------


def test_uint16_silently_truncates_beyond_its_maximum():
    """Documents the failure mode uint32 avoids: no exception, missing cells."""
    mask = _rasterize(UINT16_MAX + 2, np.uint16)

    assert mask.max() == UINT16_MAX
    assert len(np.unique(mask)) - 1 == UINT16_MAX  # two labels lost, silently


def test_label_dtype_preserves_every_label_past_uint16():
    n = UINT16_MAX + 2
    mask = _rasterize(n, LABEL_DTYPE)

    assert mask.dtype == np.uint32
    assert mask.max() == n
    assert len(np.unique(mask)) - 1 == n


# ----------------------------------------------------------------------------
# compression
# ----------------------------------------------------------------------------


def test_mask_compression_is_lossless_and_smaller(tmp_path):
    """These masks are written to and re-read from shared scratch, so raw byte
    count is the cost. Compression must not change a single label."""
    mask = _rasterize(5000, np.uint32)
    path = tmp_path / "mask.tiff"

    tifffile.imwrite(
        path, mask,
        compression=MASK_COMPRESSION, compressionargs=MASK_COMPRESSION_ARGS,
    )

    assert path.stat().st_size < mask.nbytes
    assert np.array_equal(tifffile.imread(path), mask)


# ----------------------------------------------------------------------------
# pixel convention: geometries from sopa back to the same pixels
# ----------------------------------------------------------------------------
#
# sopa vectorises a label image with pixel (row r, column c) centred on the
# point (x=c, y=r). These tests build cell outlines the same way (the union of
# the unit squares around each pixel centre) and check that rasterising them
# gives back exactly the label image they came from.


def _outlines_from_labels(labels):
    """One polygon per label: the union of its pixels' squares, centred on (c, r)."""
    from shapely.ops import unary_union

    shapes = []
    for value in np.unique(labels[labels > 0]):
        rows, cols = np.nonzero(labels == value)
        squares = [box(c - 0.5, r - 0.5, c + 0.5, r + 0.5) for r, c in zip(rows, cols)]
        shapes.append((unary_union(squares), int(value)))
    return shapes


def _example_labels():
    """Three cells in a 12 x 16 image, two of them separated by a 1-pixel gap."""
    labels = np.zeros((12, 16), dtype=np.uint32)
    labels[0:4, 0:3] = 1     # touches the top-left corner
    labels[5:10, 4:8] = 2
    labels[5:10, 9:13] = 3   # column 8 is the 1-pixel gap between cells 2 and 3
    labels[6, 13] = 3        # a one-pixel bump
    return labels


def test_rasterize_labels_reproduces_the_label_image():
    labels = _example_labels()
    height, width = labels.shape

    mask = rasterize_labels(_outlines_from_labels(labels), width, height)

    assert mask.dtype == LABEL_DTYPE
    assert np.array_equal(mask, labels)


def test_rasterize_labels_keeps_the_gap_between_neighbouring_cells():
    """all_touched=True filled the gap column, so neighbouring cells touched."""
    labels = _example_labels()
    height, width = labels.shape

    mask = rasterize_labels(_outlines_from_labels(labels), width, height)

    assert not mask[5:10, 8].any()


def test_rasterize_labels_keeps_image_orientation():
    """Row 0 of the output is y = 0 (the top of the image); no flip is needed."""
    labels = np.zeros((10, 10), dtype=np.uint32)
    labels[0, 0] = 1
    labels[9, 3] = 2

    mask = rasterize_labels(_outlines_from_labels(labels), 10, 10)

    assert mask[0, 0] == 1
    assert mask[9, 3] == 2
    assert np.count_nonzero(mask) == 2
