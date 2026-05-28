# WEHI-SODA-Hub/sp_segment: Usage

## Introduction

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated or YAML file with a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

Ensure that each row has a unique `sample` name. Even if you want to run two both Mesmer and Cellpose on the same sample, you will have to give them different sample names to run both methods without collisions.

If you only want too perform background subtraction, the minimal sample sheet is:

```csv
sample,run_backsub,tiff
sample1,true,/path/to/sample1.tiff
sample2,true,/path/to/sample2.tiff
```

### Full samplesheet

A full sample sheet is shown below:

```csv
sample,run_backsub,run_mesmer,run_cellpose,run_cellsam,tiff
sample1,true,true,false,false,/path/to/sample1.tiff
sample2,true,false,true,false,/path/to/sample2.tiff
sample3,false,false,false,true,/path/to/sample3.tiff
```

You may also prefer to use YAML for your samplesheet, either is supported:

`samplesheet.yml`:

```yaml
- sample: sample1
  run_backsub: true
  run_mesmer: true
  run_cellpose: false
  run_cellsam: false
  tiff: /path/to/sample1.tiff
- sample: sample2
  run_backsub: true
  run_mesmer: false
  run_cellpose: true
  run_cellsam: false
  tiff: /path/to/sample2.tiff
- sample: sample3
  run_backsub: false
  run_mesmer: false
  run_cellpose: false
  run_cellsam: true
  tiff: /path/to/sample3.tiff
```

| Column         | Description                                                                                                           |
| -------------- | --------------------------------------------------------------------------------------------------------------------- |
| `sample`       | Custom sample name.                                                                                                   |
| `run_backsub`  | Run background subtraction on the image                                                                               |
| `run_mesmer`   | Run Mesmer segmentation on the image (only one of `run_mesmer`, `run_cellpose`, `run_cellsam` can be true per row).   |
| `run_cellpose` | Run Cellpose segmentation on the image (only one of `run_mesmer`, `run_cellpose`, `run_cellsam` can be true per row). |
| `run_cellsam`  | Run CellSAM segmentation on the image (only one of `run_mesmer`, `run_cellpose`, `run_cellsam` can be true per row).  |
| `tiff`         | OME-TIFF for COMET or multi-channel TIFF from MIBI                                                                    |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run WEHI-SODA-Hub/sp_segment \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run WEHI-SODA-Hub/sp_segment -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: "./samplesheet.csv"
outdir: "./results/"
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Background subtraction parameters

| Parameter Name | Description                                                     |
| -------------- | --------------------------------------------------------------- |
| remove_markers | Marker channels to remove from the background subtracted image. |

### Cell processing parameters

| Parameter Name      | Description                                                                        |
| ------------------- | ---------------------------------------------------------------------------------- |
| use_whole_cell_only | Use only the whole-cell segmentation to process cells (skip nuclear segmentation). |

### Combine channel parameters

Works for both mesmer and sopa segmentation.

| Parameter Name | Description                                                |
| -------------- | ---------------------------------------------------------- |
| combine_method | Method used to combine membrane channels (max or product). |

### Mesmer parameters

The following Mesmer parameters can be set:

| Parameter Name             | Description                                                                                                                                              |
| -------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| mesmer_segmentation_level  | Segmentation level (legacy parameter).                                                                                                                   |
| mesmer_maxima_threshold    | Controls segmentation level directly in mesmer, (lower values = more cells, higher values = fewer cells).                                                |
| mesmer_interior_threshold  | Controls how conservative model is in distinguishing cell from background (higher values = more conservative interior therefore smaller nuclei & cells). |
| mesmer_maxima_smooth       | Smooths signal peaks (higher values = less irregular shapes/nuclei).                                                                                     |
| mesmer_min_nuclei_area     | Minimum area of nuclei to keep in square pixels.                                                                                                         |
| mesmer_remove_border_cells | Remove cells that touch the image border.                                                                                                                |
| mesmer_pixel_expansion     | Manual pixel expansion after segmentation.                                                                                                               |
| mesmer_padding             | Number of pixels to crop the image by on each side before segmentation.                                                                                  |

### Cellpose parameters

| Parameter Name              | Description                                                                     |
| --------------------------- | ------------------------------------------------------------------------------- |
| cellpose_diameter           | Diameter of cells in pixels for cellpose.                                       |
| cellpose_min_area           | Minimum area of cells in square pixels for cellpose.                            |
| cellpose_flow_threshold     | Flow threshold for cellpose.                                                    |
| cellpose_cellprob_threshold | Cell probability threshold for cellpose.                                        |
| cellpose_model_type         | Cellpose model to use for segmentation (e.g., nuclei, cyto, cyto2, cyto3 etc.). |
| cellpose_pretrained_model   | Path to a pre-trained Cellpose model.                                           |

### CellSAM segmentation

CellSAM is an optional segmentation backend for whole-cell and nuclear masks.
Enable it per sample with `run_cellsam: true` in the samplesheet.

For gated model downloads, set your DeepCell token as a Nextflow secret:

```bash
nextflow secrets set DEEPCELL_ACCESS_TOKEN $YOUR_TOKEN
```

If no token is provided, CellSAM uses the bundled default model.

#### CellSAM parameters

| Parameter Name                     | Description                                                                                                |
| ---------------------------------- | ---------------------------------------------------------------------------------------------------------- |
| `cellsam_bbox_threshold`           | Confidence threshold for bounding-box detections (default: `0.4`).                                         |
| `cellsam_block_size`               | Tile size in pixels used when processing large images (default: `600`).                                    |
| `cellsam_overlap`                  | Overlap in pixels between adjacent tiles (default: `250`).                                                 |
| `cellsam_iou_depth`                | Search depth (pixels) for IoU-based duplicate removal at tile boundaries (default: `250`).                 |
| `cellsam_iou_threshold`            | IoU threshold for non-maximum suppression across tiles (default: `0.5`).                                   |
| `cellsam_use_wsi`                  | Enable whole-slide-image tiling mode (default: `true`).                                                    |
| `cellsam_gauge_cell_size`          | Automatically estimate cell size from the image before segmentation (default: `false`).                    |
| `cellsam_low_contrast_enhancement` | Apply contrast enhancement before segmentation for low-contrast images (default: `false`).                 |
| `cellsam_model_path`               | Path to a custom CellSAM model checkpoint. If `null` the built-in default model is used (default: `null`). |
| `cellsam_min_area`                 | Minimum cell area in square pixels; smaller objects are discarded (default: `0`).                          |

### KRONOS embeddings

KRONOS is an optional embedding step that runs after cell measurement and
writes per-cell embeddings plus a merged GeoJSON with KRONOS features.

To enable KRONOS:

```bash
nextflow run WEHI-SODA-Hub/sp_segment \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --enable_kronos true \
   --kronos_model_path /path/to/kronos_model_dir \
   --kronos_marker_metadata /path/to/marker_metadata.csv
```

Channel names are matched case-insensitively to KRONOS marker metadata. Use
`kronos_marker_mapping` when image channel names need explicit remapping.

#### KRONOS embedding parameters

| Parameter Name              | Description                                                                                                                     |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------------- |
| `enable_kronos`             | Enable the KRONOS embedding step entirely (default: `false`).                                                                   |
| `kronos_model_path`         | **Required** path to the directory containing the pre-trained KRONOS `.pt` checkpoint.                                          |
| `kronos_marker_metadata`    | **Required** path to a CSV file mapping marker channel names to KRONOS input slots.                                             |
| `kronos_config_path`        | Path to a KRONOS YAML config file overriding default model settings (default: `null`).                                          |
| `kronos_patch_size`         | Size in pixels of the square patch extracted around each cell for embedding (default: `64`).                                    |
| `kronos_batch_size`         | Number of patches processed per inference batch (default: `32`).                                                                |
| `kronos_num_workers`        | Number of PyTorch DataLoader worker processes (default: `4`).                                                                   |
| `kronos_max_value`          | Maximum pixel intensity used for per-channel normalisation (default: `65535`).                                                  |
| `kronos_marker_mapping`     | JSON string or path mapping pipeline channel names to KRONOS marker names. Uses identity mapping when `null` (default: `null`). |
| `kronos_distance_threshold` | Maximum distance in pixels between a cell centroid and its matched nucleus; larger gaps are left unmatched (default: `5.0`).    |

### SOPA patching parameters

| Parameter Name      | Description                                                              |
| ------------------- | ------------------------------------------------------------------------ |
| technology          | Image type used for zarr conversion, only `ome_tif` is supported (COMET) |
| patch_width_pixel   | Width and height of image patch in pixels                                |
| patch_overlap_pixel | Number of pixels that image patches will overlap                         |

### Mask smoothing options

| Parameter Name     | Description                                                                                                                                                                                    |
| ------------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| smooth_masks       | Enable mask smoothing before cell measurement to reduce polygon complexity (default: `false`). Prevents StackOverflowError in QuPath's GeoJSON export for images with complex cell boundaries. |
| smooth_method      | Smoothing method: `morphological` (close+open with disk kernel, conservative) or `gaussian` (blur+threshold, stronger smoothing). Default: `morphological`.                                    |
| smooth_kernel_size | Kernel size for smoothing. For morphological: disk radius (integer). For gaussian: sigma. Larger values = more smoothing. Default: `2`.                                                        |

### Cell measurement options

| Parameter Name              | Description                                                                                                                    |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| enable_measurements         | Calculate intensity and shape measurements for cell compartments (disabling will decrease execution time)                      |
| percentiles                 | Comma-separated list of percentiles to calculate per channel. Enable measurements must be set to `true` to use this parameter. |
| pixel_size_microns          | Pixel size in microns, use 0.28 for COMET and 0.390625 for MIBI                                                                |
| estimate_cell_boundary_dist | Where no matching membrane ROI exists, expand the nucleus by this many pixels                                                  |
| dist_threshold              | Maximum centroid distance in pixels for matching a nucleus to a whole-cell ROI (default: `10.0`).                              |
| downsample_factor           | Integer downsample factor applied to image and masks before measurement, `1` = disabled (default: `1.0`).                      |
| tile_size                   | Tile size in pixels for measurement image reads (default: `2048`). Aim for 2-4x more tiles than the thread count.              |
| tile_overlap                | Tile overlap in pixels for measurement image reads (default: `200`). Set at least as large as the largest cell diameter.       |
| neighbors                   | Number of nearest neighbours for neighbourhood feature aggregation, `0` = disabled (default: `5`).                             |
| erosion_steps               | Measure intensity in 5 equal-area erosion bins from the cell/nucleus boundary inward (default: `true`).                        |
| expansion_steps             | Measure intensity in 5 equal-area expansion bins within 20 µm outward from the cell boundary (default: `true`).                |
| environment_expansion       | Measure a pericellular 20 µm environment zone around each cell (default: `true`).                                              |
| gzip_geojson                | Gzip-compress the output GeoJSON (produces `.geojson.gz`). Recommended for large whole-slide images (default: `true`).         |

### Report parameters

| Parameter Name  | Description                         |
| --------------- | ----------------------------------- |
| generate_report | Generate segmentation report for QC |

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull WEHI-SODA-Hub/sp_segment
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [WEHI-SODA-Hub/sp_segment releases page](https://github.com/WEHI-SODA-Hub/sp_segment/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility. Conda is currently not supported for the pipeline.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Not supported for this pipeline.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
