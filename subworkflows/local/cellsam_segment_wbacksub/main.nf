include { BACKGROUNDSUBTRACT } from '../backgroundsubtract/main.nf'
include { CELLSAM_SEGMENT    } from '../cellsam_segment/main.nf'

workflow CELLSAM_SEGMENT_WBACKSUB {

    take:
    ch_cellsam_wbacksub

    main:

    ch_versions = channel.empty()

    //
    // Run background subtraction
    //
    BACKGROUNDSUBTRACT(
        ch_cellsam_wbacksub.map {
            sample,
            _run_backsub,
            _run_mesmer,
            _run_cellpose,
            _run_cellsam,
            tiff,
            _nuclear_channel,
            _membrane_channels -> [
                sample,
                tiff
            ]
        }
    )
    ch_versions = ch_versions.mix(BACKGROUNDSUBTRACT.out.versions.first())


    // Replace tiff with backsub_tif
    ch_cellsam_wbacksub
        .join( BACKGROUNDSUBTRACT.out.backsub_tif )
        .map {
            sample,
            run_backsub,
            run_mesmer,
            run_cellpose,
            run_cellsam,
            _tiff,
            nuclear_channel,
            membrane_channels,
            backsub_tiff -> [
                sample,
                run_backsub,
                run_mesmer,
                run_cellpose,
                run_cellsam,
                backsub_tiff,
                nuclear_channel,
                membrane_channels
            ]
        }.set { ch_cellsam }

    CELLSAM_SEGMENT(
        ch_cellsam
    )
    ch_versions = ch_versions.mix(CELLSAM_SEGMENT.out.versions.first())

    emit:
    nuclear_segmentation_mask    = CELLSAM_SEGMENT.out.nuclear_segmentation_mask    // channel: [ val(meta), *.tiff ]
    wholecell_segmentation_mask  = CELLSAM_SEGMENT.out.wholecell_segmentation_mask  // channel: [ val(meta), *.tiff ]
    annotations                  = CELLSAM_SEGMENT.out.annotations                   // channel: [ val(meta), *.parquet ]
    kronos_embeddings            = CELLSAM_SEGMENT.out.kronos_embeddings              // channel: [ val(meta), *.csv ] OPTIONAL
    kronos_marker_report         = CELLSAM_SEGMENT.out.kronos_marker_report           // channel: [ val(meta), *.txt ] OPTIONAL
    report                       = CELLSAM_SEGMENT.out.report                        // channel: [ val(meta), *.html ]

    versions = ch_versions                                                            // channel: [ versions.yml ]
}
