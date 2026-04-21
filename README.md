<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/wehi-soda-hub-sp_segment_logo_dark.png">
    <img alt="WEHI-SODA-Hub/sp_segment" src="assets/wehi-soda-hub-sp_segment_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/WEHI-SODA-Hub/sp_segment/actions/workflows/ci.yml/badge.svg)](https://github.com/WEHI-SODA-Hub/sp_segment/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/WEHI-SODA-Hub/sp_segment/actions/workflows/linting.yml/badge.svg)](https://github.com/WEHI-SODA-Hub/sp_segment/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.17103183-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.17103183)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.04.2-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.3.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.3.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/WEHI-SODA-Hub/sp_segment)

## Introduction

**WEHI-SODA-Hub/sp_segment** is a pipeline for running cell segmentation
on COMET, MIBI, and OPAL data. For COMET, background subtraction can be performed
followed by patched cellpose segmentation, non-patched mesmer segmentation, or
CellSAM foundation model segmentation. For MIBI, mesmer or CellSAM segmentation
can be run. Whole-cell and nuclear segmentations are run separately, and then
consolidated into whole cells with nuclei with full shape and intensity
measurements per compartment. The output GeoJSON files can be viewed in QuPath.

<details>
  <summary>Click to view Mermaid diagram</summary>
  ```mermaid
    flowchart TD
    A("COMET TIFF") --> B["Extract markers"]
    B --> C["Background
            subtraction"]
    C --> D{"Segmentation
                method"} & O["Backsub TIFF"]
        N("COMET/MIBI TIFF") --> D
        D -- Cellpose (COMET only) --> S["Combine
                                         channels"]
        S --> E["sopa convert"]
        E --> F["sopa patchify"]
        F --> G["cellpose
                (nuclear)"]
        F --> H["cellpose
                (whole-cell)"]
        G --> I["sopa resolve"]
        H --> I
        I --> J["parquet to tiff"]
        J --> K["Cell measurement"]
        D -- Mesmer (COMET/MIBI) --> L["mesmer
                                       (nuclear)"]
        D -- Mesmer (COMET/MIBI) --> M["mesmer
                                       (whole-cell)"]
        L --> K
        M --> K
        K --> P("GeoJSON")
        K --> Q["segmentation
                report"]
        Q --> R("html file")
```
</details>

![sp_segment workflow](assets/sp_segment_workflow.png)

The pipeline uses the following tools:

- [Background_subtraction](https://github.com/SchapiroLabor/Background_subtraction)
  -- background subtraction tool for COMET.
- [MesmerSegmentation](https://github.com/WEHI-SODA-Hub/mesmersegmentation) -- a
  CLI for running Mesmer segmentation of MIBI and OME-XML TIFFs.
- [CellSAM](https://github.com/vanvalenlab/cellSAM) -- a foundation model for
  cell segmentation across diverse imaging modalities.
- [cellmeasurement](https://github.com/WEHI-SODA-Hub/cellmeasurement) -- a
  Python app that matches whole-cell segmentations with nuclei and calculates
  compartment measurements and intensities.
- [KRONOS](https://github.com/mahmoodlab/KRONOS) -- a foundation model for
  multiplex spatial proteomics that extracts rich embeddings for each cell.
- [sopa](https://github.com/gustaveroussy/sopa) -- we use the sopa CLI tool to
  patchify images and perform cellpose segmentation.
- [spatialVis](https://github.com/WEHI-SODA-Hub/spatialVis) -- R package for spatial
  analyses, used to generate plots for the segmentation report.

Please see the [docs for more detailed information on pipeline usage and output](docs/README.md)

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` (to test cellpose segmentation) or `-profile test_mesmer` to test mesmer segmentation before running the workflow on actual data.

If you are running this pipeline from WEHI, it has been set up to run on [Seqera Platform](https://seqera.services.biocommons.org.au/).

Usage will depend on your desired steps. See [usage docs](docs/usage.md) for more detailed information.

### Background subtraction

> [!NOTE]
> This step will only work with COMET OME-TIF files.

Prepare a sample sheet as follows:

`samplesheet.csv`:

```csv
sample,run_backsub,run_mesmer,run_cellpose,run_cellsam,tiff
sample1,true,true,false,false,/path/to/sample1.tiff
sample2,true,false,false,true,/path/to/sample2.tiff
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
  run_mesmer: true
  run_cellpose: false
  run_cellsam: false
  tiff: /path/to/sample2.tiff
```

> [!WARNING]
> Please ensure that your image name and all directories in your path do not contain spaces.

If you don't specify any segmentation algorithm to run (mesmer, cellpose, or cellsam), the pipeline will run a background subtraction only.

Now, you can run the pipeline using:

```bash
nextflow run WEHI-SODA-Hub/sp_segment \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

### Mesmer segmentation

Before running Mesmer, ensure that you have a [deepcell access token](https://users.deepcell.org/login/)
and that you have set it in your Nextflow secrets:

```bash
nextflow secrets set DEEPCELL_ACCESS_TOKEN $YOUR_TOKEN
```

If you want to run Mesmer as your segmentation algorithm, you can specify a
config file like so:

```csv
sample,run_backsub,run_mesmer,tiff,nuclear_channel,membrane_channels
sample1,true,true,/path/to/sample1.tiff,DAPI,CD45:CD8
sample2,false,true,/path/to/sample2.tiff,DAPI,CD45
```

Nuclear channels only support one entry; membrane channels may have multiple
values separated by `:` characters. If your channels have spaces in them, make
sure that you surround your channel name with quotes. For example, CD45:"HLA I".

You can also set the segmentation parameters for mesmer either via CLI
(e.g., `--combine_method prod` or in a config file pass to the workflow
via `-c`. See [usage](docs/usage.md) for a full list.

> [!NOTE]
> You cannot run multiple segmentation methods (Mesmer, Cellpose, or CellSAM) on the same sample (with the same name). If you want to run multiple methods on a sample, put it on a different line and give it a different sample name.

### Cellpose segmentation

If you want to run Cellpose as your segmentation algorithm, you can specify a
config file like so:

```csv
sample,run_backsub,run_cellpose,tiff,nuclear_channel,membrane_channels
sample1,true,true,/path/to/sample1.tiff,DAPI,CD45:CD8
sample2,false,true,/path/to/sample2.tiff,DAPI,CD45
```

As with Mesmer, nuclear channels only support one entry; membrane channels may
have multiple values separated by `:` characters. You can also set the following
parameters, either via CLI (e.g., `--combine_method prod` or in a config
file pass to the workflow via `-c`. See [usage](docs/usage.md) for a full list.

Cellpose will run in a parallelised patched workflow using sopa. To control the
patching process, you can use the `patch_width_pixel` and `patch_overlap_pixel`
parameters.

If you want to skip measurements (this may take some time for large images), you
can use set the parameter `skip_measurements` to `true`.

### KRONOS embeddings

KRONOS is a foundation model for multiplex spatial proteomics that extracts rich embeddings for each cell. These embeddings capture cellular phenotype and microenvironment context, enabling downstream analysis like clustering, classification, and spatial analysis.

To enable KRONOS embeddings:

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --skip_kronos false \
  --kronos_model_path /path/to/kronos_model \
  --kronos_marker_metadata /path/to/marker_metadata.csv \
  --kronos_merge_geojson true \
  ...
```

#### KRONOS parameters

- `--skip_kronos` (default: true): Set to `false` to enable KRONOS embedding extraction
- `--kronos_model_path` (required): Path to the KRONOS model checkpoint (.pt file)
- `--kronos_config_path` (optional): Path to KRONOS `config.json` (auto-detected from model directory if not set)
- `--kronos_marker_metadata` (required): Path to marker metadata CSV file mapping marker IDs to names
- `--kronos_merge_geojson` (default: false): Merge embeddings into the cellmeasurement GeoJSON output
- `--kronos_patch_size` (default: 64): Patch size for cell-centered crops
- `--kronos_batch_size` (default: 32): Batch size for model inference
- `--kronos_num_workers` (default: 4): Number of DataLoader workers for parallel data loading
- `--kronos_max_value` (default: 65535): Maximum intensity value for normalization
- `--kronos_marker_mapping` (optional): JSON string mapping image marker names to KRONOS marker names
- `--kronos_distance_threshold` (default: 5.0): Maximum centroid distance (pixels) used for GeoJSON merge fallback matching

#### Embeddings for filtered data with KRONOS

When `--kronos_merge_geojson` is enabled, the pipeline automatically creates a new segmentation mask directly from the GeoJSON polygons.

#### Output files

KRONOS produces the following outputs:

- `*_kronos_embeddings.csv`: CSV file with cell IDs, centroids, and 384 embedding dimensions
- `*_marker_report.txt`: Report showing which image channels were matched to KRONOS markers
- `*_kronos_merged.geojson` (if `--kronos_merge_geojson=true`): GeoJSON file with embeddings added as cell properties

The merged GeoJSON file contains all original cell measurements plus additional features (`kronos_emb_0` through `kronos_emb_#`), enabling integrated analysis of morphology, intensity, and KRONOS embeddings.

#### Marker matching

KRONOS expects specific marker names based on its training data. The pipeline automatically performs case-insensitive matching between your image channel names and the KRONOS marker metadata. For markers that don't auto-match, use `--kronos_marker_mapping`:

```bash
--kronos_marker_mapping '{"CD3e": "CD3E", "PanCK": "PANCK"}'
```

For COMET data with fluorophore suffixes in channel names, you can map them like this:

```bash
--kronos_marker_mapping '{"DAPI": "DAPI", "FOXP3_T - TRITC": "FOXP3", "CD3_T - Cy5": "CD3"}'
```

For more information about KRONOS, see the [KRONOS GitHub repository](https://github.com/mahmoodlab/KRONOS).

### CellSAM segmentation

CellSAM is a foundation model for cell segmentation that works across different
imaging modalities. To use CellSAM as your segmentation algorithm, specify a
config file like so:

```csv
sample,run_backsub,run_cellsam,tiff,nuclear_channel,membrane_channels
sample1,true,true,/path/to/sample1.tiff,DAPI,CD45:CD8
sample2,false,true,/path/to/sample2.tiff,DAPI,CD45
```

Nuclear channels only support one entry; membrane channels may have multiple
values separated by `:` characters. If your channels have spaces in them, make
sure that you surround your channel name with quotes.

CellSAM uses a tiling approach for large images and supports the following
parameters:

- `--cellsam_bbox_threshold` (default: 0.4): Confidence threshold for cell detection
- `--cellsam_block_size` (default: 600): Size of tiles for processing
- `--cellsam_overlap` (default: 56): Tile overlap for merging
- `--cellsam_iou_threshold` (default: 0.5): IOU threshold for label merging
- `--cellsam_use_wsi` (default: true): Enable tiling for large images

#### Model weights

CellSAM can automatically download the latest model weights (v1.2) from
[users.deepcell.org](https://users.deepcell.org). To use the latest weights:

1. Create an account at [users.deepcell.org](https://users.deepcell.org)
2. Generate your access token
3. Set it as a Nextflow secret:
   ```bash
   nextflow secrets set DEEPCELL_ACCESS_TOKEN $YOUR_TOKEN
   ```

If the token is not set, CellSAM will use the default bundled model weights.

> [!NOTE]
> You cannot run both Mesmer/Cellpose and CellSAM segmentation on the same sample
> (with the same name). If you want to run multiple methods on a sample, put it
> on a different line and give it a different sample name.

## Dealing with large images

You can run the pipeline with different profiles for different size images:

- `small`: for images <150GB
- `medium`: for images <300GB
- `large`: for images <600GB

> [!WARNING]
> If you are combining many membrane channels, using `prod` as the combine method
> may lead to large memory usage. In these cases, it is recommended to use `max`
> instead.

## Credits

WEHI-SODA-Hub/sp_segment was originally written by the WEHI SODA-Hub.

We thank the following people for their extensive assistance in the development of this pipeline:

- Michael McKay (@mikemcka)
- Emma Watson

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

If you use WEHI-SODA-Hub/sp_segment for your analysis, please cite it using the following doi: [10.5281/zenodo.17103183](https://doi.org/10.5281/zenodo.17103183)

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline was created using the `nf-core` template. You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
