# WEHI-SODA-Hub/sp_segment: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.4.0 - 2026-04-20

### Changed

- **cellmeasurement: Python replaces Groovy as the only implementation.** The original Groovy app
  (`https://github.com/WEHI-SODA-Hub/cellmeasurement`) is no longer used. All measurement logic
  is now in `bin/cellmeasurement.py`.

- **cellmeasurement: erosion measurement column names have changed.** The Groovy app produced
  columns named `{Channel}: {Compartment}: Eroded_{N}px: Mean/Median` where N was a fixed pixel
  depth. The Python implementation uses 5 _equal-area_ bins and names them
  `{Compartment}: ErosionBin_{N}: Mean/Median`. If you have existing Groovy-generated GeoJSON
  that you compare with new Python output, expect different column names.

- **cellmeasurement: `--erosion-steps` no longer accepts pixel-depth values.** Passing numeric
  values such as `--erosion-steps=4,7,11,14,18` (Groovy API) is not supported. Use the boolean
  flag `--erosion-steps` to enable the 5 equal-area bins.

### Added

- **KRONOS embedding output**: new `kronosembeddings` module extracts per-cell embeddings using
  the KRONOS foundation model. Produces `*_kronos_embeddings.csv` (cell IDs + 384 embedding
  dimensions), `*_marker_report.txt` (channel-to-marker match summary), and optionally
  `*_kronos_merged.geojson` (embeddings merged into the cellmeasurement GeoJSON).
  Controlled by `skip_kronos`, `kronos_model_path`, `kronos_marker_metadata`, and related
  `kronos_*` parameters. Disabled by default (`skip_kronos = true`).

- **CellSAM segmentation**: new `cellsam_segment` module adds support for the CellSAM foundation
  model as a third segmentation option alongside Mesmer and Cellpose. Supports tiled (WSI) and
  non-tiled inference via `cellsam_use_wsi`. Controlled by `cellsam_bbox_threshold`,
  `cellsam_block_size`, `cellsam_overlap`, `cellsam_iou_threshold`, and `cellsam_use_wsi`
  parameters.

- **Mask smoothing (`smooth_masks`)**: new optional `smoothmasks` module reduces polygon
  complexity. Two methods: `morphological` (disk close+open, default) and
  `shapely` (Douglas-Peucker polygon simplification). Enabled via `smooth_masks
= true`; tunable with `smooth_method` and `smooth_kernel_size`.

- `dist_threshold` pipeline parameter: controls the maximum centroid distance (pixels) for
  matching a nucleus to a whole-cell ROI in cellmeasurement (default: 10.0).
- `downsample_factor` pipeline parameter: integer downsample factor applied to image and masks
  before cellmeasurement to speed up processing on very large images (default: 1.0 = disabled).

### Removed

- `gradle_cache_dir` parameter and all references removed (leftover from the Groovy era).

## v0.3.0 - 2025-11-13

- Rename to sp_segment
- Fix issue with OME-XML processing
- Add background image preview in report
- Support for embedded report
- Add percentile calculation
- Update workflow diagram

## v0.2.0 - 2025-09-12

Add segmentation report.

## v0.1.0 - 2025-08-29

Initial release of WEHI-SODA-Hub/sp_segment, created with the [nf-core](https://nf-co.re/) template.
