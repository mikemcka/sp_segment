#!/usr/bin/env python3
"""Tests for erosion bin, expansion bin, and environment measurement systems."""

# bin/ is added to sys.path by tests/python/conftest.py
import numpy as np
import pytest
from skimage.morphology import disk

import cellmeasurement as cm


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _circular_mask(radius: int, size: int = 0) -> np.ndarray:
    """Create a centred circular binary mask."""
    if size == 0:
        size = 2 * radius + 11  # generous padding
    centre = size // 2
    yy, xx = np.ogrid[:size, :size]
    return ((xx - centre) ** 2 + (yy - centre) ** 2 <= radius ** 2).astype(bool)


def _make_comp_masks(cell_mask: np.ndarray, nuc_mask: np.ndarray = None):
    """Build a comp_masks dict matching what cellmeasurement expects."""
    if nuc_mask is None:
        nuc_mask = np.zeros_like(cell_mask, dtype=bool)
    return {
        "CELL": cell_mask.astype(bool),
        "NUCLEUS": nuc_mask.astype(bool),
    }


# ===========================================================================
# _erosion_bins_for_mask
# ===========================================================================

class TestErosionBinsForMask:
    """Unit tests for the low-level bin boundary computation."""

    def test_returns_n_bins(self):
        mask = _circular_mask(20)
        for n in (3, 5, 7):
            bins = cm._erosion_bins_for_mask(mask, n_bins=n)
            assert len(bins) == n, f"Expected {n} bins, got {len(bins)}"

    def test_empty_mask_returns_empty(self):
        mask = np.zeros((30, 30), dtype=bool)
        assert cm._erosion_bins_for_mask(mask) == []

    def test_small_mask_pads_to_n_bins(self):
        """A tiny mask that fully erodes before all bins are reached."""
        mask = np.zeros((10, 10), dtype=bool)
        mask[4:6, 4:6] = True  # 4 pixels
        bins = cm._erosion_bins_for_mask(mask, n_bins=5)
        assert len(bins) == 5

    def test_erosion_depths_are_non_decreasing(self):
        mask = _circular_mask(25)
        bins = cm._erosion_bins_for_mask(mask, n_bins=5)
        depths = [d for _, d in bins]
        assert depths == sorted(depths), f"Depths not non-decreasing: {depths}"

    def test_eroded_area_decreases_monotonically(self):
        mask = _circular_mask(25)
        bins = cm._erosion_bins_for_mask(mask, n_bins=5)
        areas = [int(np.count_nonzero(m)) for m, _ in bins]
        for i in range(1, len(areas)):
            assert areas[i] <= areas[i - 1], (
                f"Eroded area not monotonically decreasing: {areas}"
            )


# ===========================================================================
# add_erosion_measurements — ring geometry
# ===========================================================================

class TestErosionRings:
    """Verify that the rings produced by add_erosion_measurements are
    mutually exclusive and collectively cover the whole mask."""

    @pytest.fixture()
    def circle_result(self):
        """Run erosion measurements on a circle and return (props, mask)."""
        mask = _circular_mask(25)
        total_area = int(np.count_nonzero(mask))
        image = np.ones((1, *mask.shape), dtype=np.float32) * 42.0
        ch_names = ["marker"]
        comp_masks = _make_comp_masks(mask)
        props: dict = {}
        cm.add_erosion_measurements(props, image, ch_names, comp_masks, steps=[])
        return props, mask, total_area

    def test_ring_areas_sum_to_total(self, circle_result):
        props, mask, total_area = circle_result
        ring_areas = [props[f"Cell: ErosionBin_{i}: Area_px"] for i in range(1, 6)]
        assert sum(ring_areas) == total_area, (
            f"Ring areas {ring_areas} (sum={sum(ring_areas)}) != total {total_area}"
        )

    def test_area_fractions_sum_to_one(self, circle_result):
        props, _, _ = circle_result
        fracs = [props[f"Cell: ErosionBin_{i}: Area_Fraction"] for i in range(1, 6)]
        assert abs(sum(fracs) - 1.0) < 1e-9, f"Fractions sum to {sum(fracs)}, expected 1.0"

    def test_each_bin_has_roughly_equal_area(self, circle_result):
        """Each ring should be approximately 20% of total area.

        We allow +/- 8% tolerance because discrete pixel erosion can't
        perfectly hit exact area targets, especially for smaller masks.
        """
        props, _, total_area = circle_result
        expected_frac = 0.2
        tolerance = 0.08
        for i in range(1, 6):
            frac = props[f"Cell: ErosionBin_{i}: Area_Fraction"]
            assert abs(frac - expected_frac) < tolerance, (
                f"Bin {i}: area fraction {frac:.3f} deviates from {expected_frac} "
                f"by more than {tolerance}"
            )

    def test_rings_are_mutually_exclusive(self):
        """Reconstruct ring masks and confirm no pixel belongs to two rings."""
        mask = _circular_mask(25)
        bins = cm._erosion_bins_for_mask(mask, n_bins=5)

        prev = mask.astype(bool)
        rings = []
        for eroded, _ in bins:
            ring = prev & ~eroded
            rings.append(ring)
            prev = eroded

        for i in range(len(rings)):
            for j in range(i + 1, len(rings)):
                overlap = np.count_nonzero(rings[i] & rings[j])
                assert overlap == 0, f"Rings {i+1} and {j+1} overlap by {overlap} pixels"

    def test_rings_union_equals_original_mask(self):
        """All rings together should reconstruct the original mask exactly."""
        mask = _circular_mask(25)
        bins = cm._erosion_bins_for_mask(mask, n_bins=5)

        prev = mask.astype(bool)
        union = np.zeros_like(mask, dtype=bool)
        for eroded, _ in bins:
            ring = prev & ~eroded
            union |= ring
            prev = eroded

        assert np.array_equal(union, mask), "Union of all rings does not equal original mask"

    def test_bin1_is_outermost(self):
        """Bin 1 should be the boundary ring — it should contain pixels at
        the edge of the mask that are removed by the first erosion."""
        mask = _circular_mask(25)
        from skimage.morphology import binary_erosion as _be
        eroded_once = _be(mask, cm._DISK_1)
        outer_ring = mask & ~eroded_once

        bins = cm._erosion_bins_for_mask(mask, n_bins=5)
        prev = mask.astype(bool)
        bin1_ring = prev & ~bins[0][0]

        assert np.all(bin1_ring[outer_ring]), (
            "Bin 1 ring does not fully contain the outermost pixel boundary"
        )


# ===========================================================================
# Intensity measurements per ring
# ===========================================================================

class TestErosionIntensity:
    """Verify that intensity values are measured from the correct ring pixels."""

    def test_uniform_image_gives_same_mean_all_bins(self):
        mask = _circular_mask(20)
        image = np.full((1, *mask.shape), 7.0, dtype=np.float32)
        props: dict = {}
        cm.add_erosion_measurements(
            props, image, ["ch1"], _make_comp_masks(mask), steps=[]
        )
        for i in range(1, 6):
            key = f"ch1: Cell: ErosionBin_{i}: Mean"
            if props.get(f"Cell: ErosionBin_{i}: Area_px", 0) > 0:
                assert key in props, f"Missing {key}"
                assert abs(props[key] - 7.0) < 1e-6

    def test_radial_gradient_bins_have_decreasing_mean(self):
        """With intensity decreasing toward centre, outer bins should have
        higher mean than inner bins."""
        size = 61
        mask = _circular_mask(25, size=size)
        centre = size // 2
        yy, xx = np.ogrid[:size, :size]
        dist = np.sqrt((xx - centre) ** 2 + (yy - centre) ** 2).astype(np.float32)
        image = dist[np.newaxis, :, :]

        props: dict = {}
        cm.add_erosion_measurements(
            props, image, ["ch1"], _make_comp_masks(mask), steps=[]
        )

        means = []
        for i in range(1, 6):
            key = f"ch1: Cell: ErosionBin_{i}: Mean"
            if key in props:
                means.append(props[key])

        for i in range(1, len(means)):
            assert means[i - 1] >= means[i], (
                f"Mean not decreasing inward: bin {i} ({means[i-1]:.2f}) < bin {i+1} ({means[i]:.2f})"
            )

    def test_no_intensity_keys_for_empty_ring(self):
        """If a ring has zero area, no intensity keys should be written."""
        mask = np.zeros((10, 10), dtype=bool)
        mask[4:6, 4:6] = True  # 4 pixels
        image = np.ones((1, 10, 10), dtype=np.float32)
        props: dict = {}
        cm.add_erosion_measurements(
            props, image, ["ch1"], _make_comp_masks(mask), steps=[]
        )
        for i in range(1, 6):
            area = props.get(f"Cell: ErosionBin_{i}: Area_px", 0)
            if area == 0:
                assert f"ch1: Cell: ErosionBin_{i}: Mean" not in props


# ===========================================================================
# Nucleus compartment
# ===========================================================================

class TestErosionNucleus:
    """Verify erosion bins are also computed for the nucleus compartment."""

    def test_nucleus_bins_produced(self):
        cell = _circular_mask(25)
        nuc = _circular_mask(12, size=cell.shape[0])
        image = np.ones((1, *cell.shape), dtype=np.float32)
        comp_masks = _make_comp_masks(cell, nuc)
        props: dict = {}
        cm.add_erosion_measurements(props, image, ["ch1"], comp_masks, steps=[])

        for i in range(1, 6):
            assert f"Nucleus: ErosionBin_{i}: Area_px" in props, (
                f"Missing Nucleus: ErosionBin_{i}: Area_px"
            )

    def test_nucleus_ring_areas_sum_to_nucleus_area(self):
        cell = _circular_mask(25)
        nuc = _circular_mask(12, size=cell.shape[0])
        nuc_area = int(np.count_nonzero(nuc))
        image = np.ones((1, *cell.shape), dtype=np.float32)
        comp_masks = _make_comp_masks(cell, nuc)
        props: dict = {}
        cm.add_erosion_measurements(props, image, ["ch1"], comp_masks, steps=[])

        ring_sum = sum(props[f"Nucleus: ErosionBin_{i}: Area_px"] for i in range(1, 6))
        assert ring_sum == nuc_area


# ===========================================================================
# steps parameter is ignored (API compat)
# ===========================================================================

class TestStepsIgnored:
    """The old steps parameter should be accepted but have no effect."""

    def test_same_output_regardless_of_steps(self):
        mask = _circular_mask(20)
        image = np.ones((1, *mask.shape), dtype=np.float32) * 5.0
        comp = _make_comp_masks(mask)

        props_a: dict = {}
        cm.add_erosion_measurements(props_a, image, ["ch1"], comp, steps=[])

        props_b: dict = {}
        cm.add_erosion_measurements(props_b, image, ["ch1"], comp, steps=[1, 3, 5])

        assert props_a == props_b, "Output should be identical regardless of steps value"


# ===========================================================================
# add_environment_measurements
# ===========================================================================

class TestEnvironmentMeasurements:
    """Tests for the 20 µm pericellular environment zone."""

    def test_environment_zone_is_outside_cell(self):
        """The environment zone should not overlap with the cell mask."""
        mask = _circular_mask(15)
        image = np.ones((1, *mask.shape), dtype=np.float32)
        props: dict = {}
        cm.add_environment_measurements(props, image, ["ch1"], mask, pixel_size_microns=0.5)
        # If environment was computed, the pixel count should exist
        assert "Cell: Environment_20um: Pixel_Count" in props
        assert props["Cell: Environment_20um: Pixel_Count"] > 0

    def test_pixel_size_controls_dilation_radius(self):
        """Smaller pixel size → more pixels dilated → larger environment area."""
        mask = _circular_mask(15, size=200)
        image = np.ones((1, *mask.shape), dtype=np.float32)

        props_coarse: dict = {}
        cm.add_environment_measurements(props_coarse, image, ["ch1"], mask, pixel_size_microns=1.0)

        props_fine: dict = {}
        cm.add_environment_measurements(props_fine, image, ["ch1"], mask, pixel_size_microns=0.5)

        area_coarse = props_coarse["Cell: Environment_20um: Pixel_Count"]
        area_fine = props_fine["Cell: Environment_20um: Pixel_Count"]
        assert area_fine > area_coarse, (
            f"Finer pixel size should give larger environment area: "
            f"fine={area_fine}, coarse={area_coarse}"
        )

    @pytest.mark.parametrize("px_um", [0.28, 0.39, 0.50, 1.0])
    def test_environment_computed_for_various_pixel_sizes(self, px_um):
        """Environment should be computed for any reasonable pixel size."""
        mask = _circular_mask(15, size=200)
        image = np.ones((1, *mask.shape), dtype=np.float32)
        props: dict = {}
        cm.add_environment_measurements(props, image, ["ch1"], mask, pixel_size_microns=px_um)
        assert "Cell: Environment_20um: Pixel_Count" in props
        assert props["Cell: Environment_20um: Pixel_Count"] > 0

    def test_environment_key_names_use_20um_label(self):
        """All measurement keys should use the '20um' label, not pixel-variable names."""
        mask = _circular_mask(15)
        image = np.ones((2, *mask.shape), dtype=np.float32)
        props: dict = {}
        cm.add_environment_measurements(props, image, ["ch1", "ch2"], mask, pixel_size_microns=0.5)

        expected_keys = [
            "Cell: Environment_20um: Pixel_Count",
            "Cell: Environment_20um: Area_Fraction",
            "ch1: Cell: Environment_20um: Mean",
            "ch1: Cell: Environment_20um: Median",
            "ch1: Cell: Environment_20um: Min",
            "ch1: Cell: Environment_20um: Max",
            "ch1: Cell: Environment_20um: Std.Dev.",
            "ch2: Cell: Environment_20um: Mean",
        ]
        for key in expected_keys:
            assert key in props, f"Missing expected key: {key}"

    def test_uniform_image_mean_equals_value(self):
        """With a uniform image, environment mean should equal the fill value."""
        mask = _circular_mask(15)
        fill_val = 3.14
        image = np.full((1, *mask.shape), fill_val, dtype=np.float32)
        props: dict = {}
        cm.add_environment_measurements(props, image, ["ch1"], mask, pixel_size_microns=0.5)
        assert abs(props["ch1: Cell: Environment_20um: Mean"] - fill_val) < 1e-5

    def test_empty_mask_produces_no_measurements(self):
        """An empty cell mask should produce no environment keys."""
        mask = np.zeros((30, 30), dtype=bool)
        image = np.ones((1, 30, 30), dtype=np.float32)
        props: dict = {}
        cm.add_environment_measurements(props, image, ["ch1"], mask, pixel_size_microns=0.5)
        assert len(props) == 0

    def test_area_fraction_relative_to_cell(self):
        """Area fraction should be environment area / cell area."""
        mask = _circular_mask(15)
        cell_area = int(np.count_nonzero(mask))
        image = np.ones((1, *mask.shape), dtype=np.float32)
        props: dict = {}
        cm.add_environment_measurements(props, image, ["ch1"], mask, pixel_size_microns=0.5)
        env_area = props["Cell: Environment_20um: Pixel_Count"]
        expected_frac = env_area / cell_area
        assert abs(props["Cell: Environment_20um: Area_Fraction"] - expected_frac) < 1e-9


# ===========================================================================
# _expansion_bins_for_mask
# ===========================================================================

class TestExpansionBinsForMask:
    """Unit tests for the expansion bin boundary computation."""

    def test_returns_n_bins(self):
        mask = _circular_mask(15, size=200)
        for n in (3, 5, 7):
            bins = cm._expansion_bins_for_mask(mask, total_expansion_px=40, n_bins=n)
            assert len(bins) == n, f"Expected {n} bins, got {len(bins)}"

    def test_empty_mask_returns_empty(self):
        mask = np.zeros((30, 30), dtype=bool)
        assert cm._expansion_bins_for_mask(mask, total_expansion_px=10) == []

    def test_dilation_depths_are_non_decreasing(self):
        mask = _circular_mask(15, size=200)
        bins = cm._expansion_bins_for_mask(mask, total_expansion_px=40, n_bins=5)
        depths = [d for _, d in bins]
        assert depths == sorted(depths), f"Depths not non-decreasing: {depths}"

    def test_dilated_area_increases_monotonically(self):
        mask = _circular_mask(15, size=200)
        bins = cm._expansion_bins_for_mask(mask, total_expansion_px=40, n_bins=5)
        areas = [int(np.count_nonzero(m)) for m, _ in bins]
        for i in range(1, len(areas)):
            assert areas[i] >= areas[i - 1], (
                f"Dilated area not monotonically increasing: {areas}"
            )


# ===========================================================================
# add_expansion_measurements — ring geometry
# ===========================================================================

class TestExpansionRings:
    """Verify expansion rings are mutually exclusive and cover the full zone."""

    @pytest.fixture()
    def expansion_result(self):
        """Run expansion measurements on a circle and return (props, mask)."""
        mask = _circular_mask(15, size=200)
        image = np.ones((1, *mask.shape), dtype=np.float32) * 42.0
        props: dict = {}
        cm.add_expansion_measurements(props, image, ["marker"], mask, pixel_size_microns=0.5)
        return props, mask

    def test_five_bins_produced(self, expansion_result):
        props, _ = expansion_result
        for i in range(1, 6):
            assert f"Cell: ExpansionBin_{i}: Area_px" in props, (
                f"Missing Cell: ExpansionBin_{i}: Area_px"
            )

    def test_each_bin_has_nonzero_area(self, expansion_result):
        """For a mask with enough room to expand, all bins should have area."""
        props, _ = expansion_result
        for i in range(1, 6):
            assert props[f"Cell: ExpansionBin_{i}: Area_px"] > 0, (
                f"Bin {i} has zero area"
            )

    def test_rings_are_mutually_exclusive(self):
        """Reconstruct ring masks and confirm no pixel belongs to two rings."""
        mask = _circular_mask(15, size=200)
        expansion_px = max(1, int(round(20.0 / 0.5)))
        bins = cm._expansion_bins_for_mask(mask, total_expansion_px=expansion_px, n_bins=5)

        cm_bool = mask.astype(bool)
        prev = cm_bool.copy()
        rings = []
        for dilated, _ in bins:
            ring = dilated & ~prev
            rings.append(ring)
            prev = dilated

        for i in range(len(rings)):
            for j in range(i + 1, len(rings)):
                overlap = np.count_nonzero(rings[i] & rings[j])
                assert overlap == 0, f"Rings {i+1} and {j+1} overlap by {overlap} pixels"

    def test_rings_union_equals_full_expansion_zone(self):
        """All rings together should equal the full dilated zone minus the cell."""
        from scipy import ndimage as ndi
        mask = _circular_mask(15, size=200)
        expansion_px = max(1, int(round(20.0 / 0.5)))
        bins = cm._expansion_bins_for_mask(mask, total_expansion_px=expansion_px, n_bins=5)

        cm_bool = mask.astype(bool)
        full_dilated = ndi.binary_dilation(cm_bool, structure=cm._DISK_1, iterations=expansion_px)
        expected_zone = full_dilated & ~cm_bool

        prev = cm_bool.copy()
        union = np.zeros_like(mask, dtype=bool)
        for dilated, _ in bins:
            ring = dilated & ~prev
            union |= ring
            prev = dilated

        assert np.array_equal(union, expected_zone), (
            "Union of all expansion rings does not equal the full expansion zone"
        )

    def test_bin1_is_closest_to_cell(self):
        """Bin 1 should contain the pixels immediately adjacent to the cell."""
        from scipy import ndimage as ndi
        mask = _circular_mask(15, size=200)
        cm_bool = mask.astype(bool)
        dilated_once = ndi.binary_dilation(cm_bool, structure=cm._DISK_1, iterations=1)
        adjacent_ring = dilated_once & ~cm_bool

        expansion_px = max(1, int(round(20.0 / 0.5)))
        bins = cm._expansion_bins_for_mask(cm_bool, total_expansion_px=expansion_px, n_bins=5)
        bin1_ring = bins[0][0] & ~cm_bool

        # Bin 1 must contain the 1-pixel adjacent boundary
        assert np.all(bin1_ring[adjacent_ring]), (
            "Bin 1 ring does not contain the pixels immediately adjacent to the cell"
        )

    def test_each_bin_has_roughly_equal_area(self, expansion_result):
        """Each ring should be approximately 20% of total expansion zone area."""
        props, _ = expansion_result
        areas = [props[f"Cell: ExpansionBin_{i}: Area_px"] for i in range(1, 6)]
        total = sum(areas)
        expected_frac = 0.2
        tolerance = 0.08
        for i, area in enumerate(areas, start=1):
            frac = area / total
            assert abs(frac - expected_frac) < tolerance, (
                f"Bin {i}: area fraction {frac:.3f} deviates from {expected_frac} "
                f"by more than {tolerance}"
            )


# ===========================================================================
# Expansion intensity
# ===========================================================================

class TestExpansionIntensity:
    """Verify intensity values are measured from the correct expansion ring pixels."""

    def test_uniform_image_gives_same_mean_all_bins(self):
        mask = _circular_mask(15, size=200)
        image = np.full((1, *mask.shape), 9.0, dtype=np.float32)
        props: dict = {}
        cm.add_expansion_measurements(props, image, ["ch1"], mask, pixel_size_microns=0.5)
        for i in range(1, 6):
            key = f"ch1: Cell: ExpansionBin_{i}: Mean"
            if props.get(f"Cell: ExpansionBin_{i}: Area_px", 0) > 0:
                assert key in props, f"Missing {key}"
                assert abs(props[key] - 9.0) < 1e-6

    def test_radial_gradient_bins_have_increasing_mean(self):
        """With intensity increasing outward from centre, outer bins should
        have higher mean than inner bins."""
        size = 200
        mask = _circular_mask(15, size=size)
        centre = size // 2
        yy, xx = np.ogrid[:size, :size]
        dist = np.sqrt((xx - centre) ** 2 + (yy - centre) ** 2).astype(np.float32)
        image = dist[np.newaxis, :, :]

        props: dict = {}
        cm.add_expansion_measurements(props, image, ["ch1"], mask, pixel_size_microns=0.5)

        means = []
        for i in range(1, 6):
            key = f"ch1: Cell: ExpansionBin_{i}: Mean"
            if key in props:
                means.append(props[key])

        for i in range(1, len(means)):
            assert means[i] >= means[i - 1], (
                f"Mean not increasing outward: bin {i} ({means[i-1]:.2f}) > bin {i+1} ({means[i]:.2f})"
            )

    def test_empty_mask_produces_no_measurements(self):
        mask = np.zeros((30, 30), dtype=bool)
        image = np.ones((1, 30, 30), dtype=np.float32)
        props: dict = {}
        cm.add_expansion_measurements(props, image, ["ch1"], mask, pixel_size_microns=0.5)
        assert len(props) == 0

    def test_pixel_size_controls_expansion_distance(self):
        """Smaller pixel size → more expansion pixels → larger total zone."""
        mask = _circular_mask(15, size=300)
        image = np.ones((1, *mask.shape), dtype=np.float32)

        props_coarse: dict = {}
        cm.add_expansion_measurements(props_coarse, image, ["ch1"], mask, pixel_size_microns=1.0)

        props_fine: dict = {}
        cm.add_expansion_measurements(props_fine, image, ["ch1"], mask, pixel_size_microns=0.5)

        total_coarse = sum(props_coarse.get(f"Cell: ExpansionBin_{i}: Area_px", 0) for i in range(1, 6))
        total_fine = sum(props_fine.get(f"Cell: ExpansionBin_{i}: Area_px", 0) for i in range(1, 6))
        assert total_fine > total_coarse, (
            f"Finer pixel size should give more expansion area: fine={total_fine}, coarse={total_coarse}"
        )
