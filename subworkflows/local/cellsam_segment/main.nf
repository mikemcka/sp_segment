include { CELLSAMSEGMENT as CELLSAMWC    } from '../../../modules/local/cellsamsegment/main.nf'
include { CELLSAMSEGMENT as CELLSAMNUC   } from '../../../modules/local/cellsamsegment/main.nf'
include { CELLMEASUREMENT                } from '../../../modules/local/cellmeasurement/main.nf'
include { SMOOTHMASKS as SMOOTHMASKS_NUC } from '../../../modules/local/smoothmasks/main.nf'
include { SMOOTHMASKS as SMOOTHMASKS_WC  } from '../../../modules/local/smoothmasks/main.nf'
include { KRONOSEMBEDDINGS               } from '../../../modules/local/kronosembeddings/main.nf'
include { COMBINECHANNELS                } from '../../../modules/local/combinechannels/main.nf'
include { SEGMENTATIONREPORT             } from '../../../modules/local/segmentationreport/main.nf'

workflow CELLSAM_SEGMENT {

    take:
    ch_cellsam_segment // channel: [ (sample, run_backsub, run_mesmer, run_cellpose, run_cellsam, tiff, nuclear_channel, membrane_channels) ]

    main:

    ch_versions = channel.empty()

    ch_cellsam_segment.map {
        sample,
        _run_backsub,
        _run_mesmer,
        _run_cellpose,
        _run_cellsam,
        tiff,
        nuclear_channel,
        membrane_channels -> [
            sample,
            tiff,
            nuclear_channel,
            membrane_channels
        ]
    }.set { ch_cellsam }


    //
    // Run CELLSAMSEGMENT module for whole-cell segmentation
    //
    CELLSAMWC(
        ch_cellsam,
        "whole-cell"
    )
    ch_versions = ch_versions.mix(CELLSAMWC.out.versions.first())


    //
    // Run CELLSAMSEGMENT module for nuclear segmentation (skipped if use_whole_cell_only is true)
    //
    if (!params.use_whole_cell_only) {
        CELLSAMNUC(
            ch_cellsam,
            "nuclear"
        )
        ch_versions = ch_versions.mix(CELLSAMNUC.out.versions.first())
    }

    // Create channel for CELLMEASUREMENT input adding the segmentation masks
    if (params.use_whole_cell_only) {
        // In whole-cell-only mode, provide whole-cell mask as a placeholder for the
        // nuclear input channel; CELLMEASUREMENT omits --nuclear-mask in this mode.
        ch_cellsam_segment
            .join(CELLSAMWC.out.segmentation_mask)
            .map {
                sample,
                _run_backsub,
                _run_mesmer,
                _run_cellpose,
                _run_cellsam,
                tiff,
                _nuclear_channel,
                _membrane_channels,
                whole_cell_mask -> [
                    sample,
                    tiff,
                    whole_cell_mask,  // placeholder for nuclear mask
                    whole_cell_mask
                ]
            }.set { ch_cellmeasurement }
    } else {
        ch_cellsam_segment
            .join(CELLSAMNUC.out.segmentation_mask)
            .join(CELLSAMWC.out.segmentation_mask)
            .map {
                sample,
                _run_backsub,
                _run_mesmer,
                _run_cellpose,
                _run_cellsam,
                tiff,
                _nuclear_channel,
                _membrane_channels,
                nuclear_mask,
                whole_cell_mask -> [
                    sample,
                    tiff,
                    nuclear_mask,
                    whole_cell_mask
                ]
            }.set { ch_cellmeasurement }
    }

    //
    // Optional mask smoothing to reduce polygon complexity
    //
    if (params.smooth_masks) {
        if (!params.use_whole_cell_only) {
            SMOOTHMASKS_NUC(
                ch_cellmeasurement.map { sample, _tiff, nuclear_mask, _whole_cell_mask -> [sample, nuclear_mask] }
            )
        }
        SMOOTHMASKS_WC(
            ch_cellmeasurement.map { sample, _tiff, _nuclear_mask, whole_cell_mask -> [sample, whole_cell_mask] }
        )
        if (!params.use_whole_cell_only) {
            ch_cellmeasurement
                .map { sample, tiff, _nuclear_mask, _whole_cell_mask -> [sample, tiff] }
                .join(SMOOTHMASKS_NUC.out.smoothed_mask)
                .join(SMOOTHMASKS_WC.out.smoothed_mask)
                .set { ch_cellmeasurement }
            ch_versions = ch_versions.mix(SMOOTHMASKS_NUC.out.versions.first())
        } else {
            ch_cellmeasurement
                .map { sample, tiff, _nuclear_mask, _whole_cell_mask -> [sample, tiff] }
                .join(SMOOTHMASKS_WC.out.smoothed_mask)
                .map { sample, tiff, smoothed_wc -> [sample, tiff, smoothed_wc, smoothed_wc] }
                .set { ch_cellmeasurement }
        }
        ch_versions = ch_versions.mix(SMOOTHMASKS_WC.out.versions.first())
    }

    //
    // Run CELLMEASUREMENT module on the whole-cell and nuclear segmentation masks
    //
    CELLMEASUREMENT(
        ch_cellmeasurement
    )
    ch_versions = ch_versions.mix(CELLMEASUREMENT.out.versions.first())

    ch_annotations = CELLMEASUREMENT.out.annotations

    //
    // Optional KRONOS embedding extraction
    //
    ch_kronos_embeddings = channel.empty()
    ch_kronos_marker_report = channel.empty()
    if (params.enable_kronos) {

        // Create channel for KRONOS input: tiff + whole-cell mask + geojson
        ch_cellsam_segment
            .join(CELLSAMWC.out.segmentation_mask)
            .map {
                sample,
                _run_backsub,
                _run_mesmer,
                _run_cellpose,
                _run_cellsam,
                tiff,
                _nuclear_channel,
                _membrane_channels,
                whole_cell_mask -> [
                    sample,
                    tiff,
                    whole_cell_mask
                ]
            }.set { ch_kronos_input }

        KRONOSEMBEDDINGS(
            ch_kronos_input,
            file(params.kronos_model_path),
            file(params.kronos_marker_metadata),
            CELLMEASUREMENT.out.annotations
        )
        ch_versions = ch_versions.mix(KRONOSEMBEDDINGS.out.versions.first())
        ch_kronos_embeddings = KRONOSEMBEDDINGS.out.embeddings
        ch_kronos_marker_report = KRONOSEMBEDDINGS.out.marker_report
        ch_annotations = KRONOSEMBEDDINGS.out.merged_geojson
    }

    // Optional SEGMENTATIONREPORT module
    ch_report = channel.empty()
    if (params.generate_report) {

        //
        // Combine channels for report background image
        //
        COMBINECHANNELS(
            ch_cellsam
        )
        ch_versions = ch_versions.mix(COMBINECHANNELS.out.versions.first())

        ch_cellsam_segment
            .join(ch_annotations)
            .join(COMBINECHANNELS.out.combined_tiff, by: 0)
            .map {
                sample,
                _run_backsub,
                run_mesmer,
                run_cellpose,
                run_cellsam,
                _tiff,
                nuclear_channel,
                membrane_channels,
                annotations,
                combined_tiff -> [
                    sample,
                    annotations,
                    run_mesmer,
                    (run_cellpose || run_cellsam),
                    nuclear_channel,
                    membrane_channels,
                    combined_tiff
                ]
            }.set { ch_segmentation_report }

        //
        // Generate segmentation report
        //
        SEGMENTATIONREPORT(
            ch_segmentation_report
        )
        ch_versions = ch_versions.mix(SEGMENTATIONREPORT.out.versions.first())
        ch_report = SEGMENTATIONREPORT.out.report
    }

    emit:
    nuclear_segmentation_mask    = params.use_whole_cell_only ? channel.empty() : CELLSAMNUC.out.segmentation_mask  // channel: [ val(meta), *.tiff ]
    wholecell_segmentation_mask  = CELLSAMWC.out.segmentation_mask        // channel: [ val(meta), *.tiff ]
    annotations                  = ch_annotations                         // channel: [ val(meta), *.geojson ]
    kronos_embeddings            = ch_kronos_embeddings                   // channel: [ val(meta), *.csv ] OPTIONAL
    kronos_marker_report         = ch_kronos_marker_report                // channel: [ val(meta), *.txt ] OPTIONAL
    report                       = ch_report                              // channel: [ val(meta), *.html ]

    versions = ch_versions                                                // channel: [ versions.yml ]
}
