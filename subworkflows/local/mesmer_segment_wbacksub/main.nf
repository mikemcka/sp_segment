include { BACKGROUNDSUBTRACT } from '../backgroundsubtract/main.nf'
include { MESMER_SEGMENT     } from '../mesmer_segment/main.nf'

workflow MESMER_SEGMENT_WBACKSUB {

    take:
    ch_mesmer_wbacksub

    main:

    ch_versions = channel.empty()

    //
    // Run background subtraction
    //
    BACKGROUNDSUBTRACT(
        ch_mesmer_wbacksub.map {
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
    ch_mesmer_wbacksub
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
        }.set { ch_mesmer }

    MESMER_SEGMENT(
        ch_mesmer
    )
    ch_versions = ch_versions.mix(MESMER_SEGMENT.out.versions.first())

    emit:
    annotations      = MESMER_SEGMENT.out.annotations       // channel: [ val(meta), *.geojson ]
    whole_cell_tif   = MESMER_SEGMENT.out.whole_cell_tif    // channel: [ val(meta), *.tiff ]
    nuclear_tif      = MESMER_SEGMENT.out.nuclear_tif       // channel: [ val(meta), *.tiff ]
    kronos_embeddings     = MESMER_SEGMENT.out.kronos_embeddings   // channel: [ val(meta), *.csv ] OPTIONAL
    kronos_marker_report  = MESMER_SEGMENT.out.kronos_marker_report // channel: [ val(meta), *.txt ] OPTIONAL

    versions = ch_versions                                  // channel: [ versions.yml ]
}
