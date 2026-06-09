include { COMBINECHANNELS                                    } from '../../../modules/local/combinechannels/main.nf'
include { SOPA_SEGMENT_COMPARTMENT as SOPA_SEGMENT_NUCLEAR   } from '../sopa_segment_compartment/main.nf'
include { SOPA_SEGMENT_COMPARTMENT as SOPA_SEGMENT_WHOLECELL } from '../sopa_segment_compartment/main.nf'
include { CELLMEASUREMENT                                    } from '../../../modules/local/cellmeasurement/main.nf'
include { SMOOTHMASKS as SMOOTHMASKS_NUC                     } from '../../../modules/local/smoothmasks/main.nf'
include { SMOOTHMASKS as SMOOTHMASKS_WC                      } from '../../../modules/local/smoothmasks/main.nf'
include { KRONOSEMBEDDINGS                                    } from '../../../modules/local/kronosembeddings/main.nf'
include { SEGMENTATIONREPORT                                 } from '../../../modules/local/segmentationreport/main.nf'

workflow SOPA_SEGMENT {

    take:
    ch_sopa // channel: [ (meta, tiff, nuclear_channel, membrane_channels) ]

    main:

    ch_versions = channel.empty()

    //
    // Combine membrane channels into a single channel
    //
    COMBINECHANNELS(
        ch_sopa
    )
    ch_versions = ch_versions.mix(COMBINECHANNELS.out.versions.first())

    // Replace tiff with combined_tiff
    // If there are multiple membrane channels, rename to 'combined_membrane'
    COMBINECHANNELS.out.combined_tiff
        .join( ch_sopa, by: 0 )
        .map { meta, combined_tiff, _tiff, nuclear_channel, membrane_channels ->
            def membrane_name = membrane_channels.split(':').size() == 1 ?
                membrane_channels : 'combined_membrane'
            [ meta, combined_tiff, nuclear_channel, membrane_name ]
        }.set { ch_combined }

    //
    // Run segmentation for nuclear compartment (skipped if use_whole_cell_only is true)
    //
    if (!params.use_whole_cell_only) {
        SOPA_SEGMENT_NUCLEAR(
            ch_combined.map {
                meta,
                tiff,
                nuclear_channel,
                _membrane_channels -> [
                    meta,
                    tiff,
                    nuclear_channel,
                    '' // no membrane channels for nuclear segmentation
                ]
            },
            'nuclear'
        )
        ch_versions = ch_versions.mix(SOPA_SEGMENT_NUCLEAR.out.versions.first())
    }

    //
    // Run segmentation for whole-cell compartment
    //
    SOPA_SEGMENT_WHOLECELL(
        ch_combined,
        'whole-cell'
    )
    ch_versions = ch_versions.mix(SOPA_SEGMENT_WHOLECELL.out.versions.first())

    //
    // Create a channel for cell measurement
    //
    if (params.use_whole_cell_only) {
        // In whole-cell-only mode, provide whole-cell mask as a placeholder for the
        // nuclear input channel; CELLMEASUREMENT omits --nuclear-mask in this mode.
        SOPA_SEGMENT_WHOLECELL.out.tiff
            .join(ch_sopa, by: 0)
            .map {
                meta,
                wholecell_tiff,
                tiff,
                _nuc_chan,
                _mem_chans -> [
                    meta,
                    tiff,
                    wholecell_tiff,  // placeholder for nuclear mask
                    wholecell_tiff
                ]
            }.set { ch_cellmeasurement }
    } else {
        SOPA_SEGMENT_NUCLEAR.out.tiff
            .join(SOPA_SEGMENT_WHOLECELL.out.tiff, by: 0)
            .join(ch_sopa, by: 0)
            .map {
                meta,
                nuclear_tiff,
                wholecell_tiff,
                tiff,
                _nuc_chan,
                _mem_chans -> [
                    meta,
                    tiff,
                    nuclear_tiff,
                    wholecell_tiff
                ]
            }.set { ch_cellmeasurement }
    }

    //
    // Optional mask smoothing to reduce polygon complexity
    //
    if (params.smooth_masks) {
        if (!params.use_whole_cell_only) {
            SMOOTHMASKS_NUC(
                ch_cellmeasurement.map {
                    meta, _tiff, nuclear_tiff, _wholecell_tiff -> [meta, nuclear_tiff]
                }
            )
        }
        SMOOTHMASKS_WC(
            ch_cellmeasurement.map {
                meta, _tiff, _nuclear_tiff, wholecell_tiff -> [meta, wholecell_tiff]
            }
        )
        if (!params.use_whole_cell_only) {
            ch_cellmeasurement
                .map { meta, tiff, _nuclear_tiff, _wholecell_tiff -> [meta, tiff] }
                .join(SMOOTHMASKS_NUC.out.smoothed_mask)
                .join(SMOOTHMASKS_WC.out.smoothed_mask)
                .set { ch_cellmeasurement }
            ch_versions = ch_versions.mix(SMOOTHMASKS_NUC.out.versions.first())
        } else {
            ch_cellmeasurement
                .map { meta, tiff, _nuclear_tiff, _wholecell_tiff -> [meta, tiff] }
                .join(SMOOTHMASKS_WC.out.smoothed_mask)
                .map { meta, tiff, smoothed_wc -> [meta, tiff, smoothed_wc, smoothed_wc] }
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

        // Create channel for KRONOS input: original tiff + whole-cell mask + geojson
        ch_sopa
            .join(SOPA_SEGMENT_WHOLECELL.out.tiff)
            .map {
                meta,
                tiff,
                _nuclear_channel,
                _membrane_channels,
                whole_cell_mask -> [
                    meta,
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
        ch_combined
            .join(ch_annotations)
            .map {
                sample,
                combined_tiff,
                nuclear_channel,
                membrane_channels,
                annotations -> [
                    sample,
                    annotations,
                    false, // run_mesmer
                    true,  // run_cellpose
                    false, // run_cellsam
                    nuclear_channel,
                    membrane_channels,
                    combined_tiff
                ]
            }.set { ch_segmentationreport }

        //
        // Run SEGMENTATIONREPORT module to generate a report of the segmentation results
        //
        SEGMENTATIONREPORT(
            ch_segmentationreport
        )
        ch_versions = ch_versions.mix(SEGMENTATIONREPORT.out.versions.first())

        ch_report = SEGMENTATIONREPORT.out.report  // channel: [ val(meta), *.html ]
    }

    emit:
    annotations          = ch_annotations                    // channel: [ val(meta), *.geojson ]
    kronos_embeddings    = ch_kronos_embeddings              // channel: [ val(meta), *.csv ] OPTIONAL
    kronos_marker_report = ch_kronos_marker_report           // channel: [ val(meta), *.txt ] OPTIONAL
    report               = ch_report                         // channel: [ val(meta), *.html ] OPTIONAL

    versions = ch_versions                                   // channel: [ versions.yml ]
}
