/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { paramsSummaryMap          } from 'plugin/nf-schema'
include { BACKGROUNDSUBTRACT        } from '../subworkflows/local/backgroundsubtract'
include { MESMER_SEGMENT_WBACKSUB   } from '../subworkflows/local/mesmer_segment_wbacksub'
include { MESMER_SEGMENT            } from '../subworkflows/local/mesmer_segment'
include { CELLSAM_SEGMENT_WBACKSUB  } from '../subworkflows/local/cellsam_segment_wbacksub'
include { CELLSAM_SEGMENT           } from '../subworkflows/local/cellsam_segment'
include { SOPA_SEGMENT              } from '../subworkflows/local/sopa_segment'
include { SOPA_SEGMENT_WBACKSUB     } from '../subworkflows/local/sopa_segment_wbacksub'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_sp_segment_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SP_SEGMENT {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()

    //
    // Construct channel for background subtraction/segmentation workflow for MESMER
    //
    ch_samplesheet.branch { it ->
        backsub_only: it[1] == true &&  // run_backsub true
                        it[2] == false && // run_mesmer false
                        it[3] == false && // run_cellpose false
                        it[4] == false    // run_cellsam false
        backsub_mesmer: it[1] == true && it[2] == true // run_backsub true, run_mesmer true
        mesmer_only: it[1] == false && it[2] == true // run_backsub false, run_mesmer true
    }.set { ch_mesmer }

    //
    // Run the BACKGROUNDSUBTRACT subworkflow for samples that ONLY require
    // background subtraction (no segmentation)
    //
    BACKGROUNDSUBTRACT(
        ch_mesmer.backsub_only.map {
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

    //
    // Run the MESMER_SEGMENT_WBACKSUB subworkflow for samples that require
    // background subtraction and mesmer segmentation
    //
    MESMER_SEGMENT_WBACKSUB(
        ch_mesmer.backsub_mesmer
    )

    //
    // Run MESMER_SEGMENT subworkflow for samples that ONLY require mesmer segmentation
    //
    MESMER_SEGMENT(
        ch_mesmer.mesmer_only
    )

    //
    // Construct channel for only CELLPOSE subworkflow
    //
    ch_samplesheet.filter {
        it[3] == true // run_cellpose true for sample
    }.map {
        sample,
        run_backsub,
        _run_mesmer,
        _run_cellpose,
        _run_cellsam,
        tiff,
        nuclear_channel,
        membrane_channels -> [
            sample,
            run_backsub,
            tiff,
            nuclear_channel,
            membrane_channels
        ]
    }.branch { it ->
        with_backsub: it[1] == true// run_backsub true
        no_backsub: it[1] == false // run_backsub false
    }.set { ch_cellpose_samplesheet }

    //
    // Run CELLPOSE subworkflow for samples that require background subtraction
    //
    SOPA_SEGMENT_WBACKSUB(
        ch_cellpose_samplesheet.with_backsub.map { sample,
            _run_backsub,
            tiff,
            nuclear_channel,
            membrane_channels ->
            [ sample, tiff, nuclear_channel, membrane_channels ]
        }
    )

    //
    // Run CELLPOSE subworkflow for samples that ONLY require cellpose segmentation
    //
    SOPA_SEGMENT(
        ch_cellpose_samplesheet.no_backsub.map { sample,
            _run_backsub,
            tiff,
            nuclear_channel,
            membrane_channels ->
            [ sample, tiff, nuclear_channel, membrane_channels ]
        }
    )

    //
    // Construct channel for CellSAM segmentation workflow
    //
    ch_samplesheet.branch { it ->
        backsub_cellsam: it[1] == true && it[4] == true // run_backsub true, run_cellsam true
        cellsam_only: it[1] == false && it[4] == true   // run_backsub false, run_cellsam true
    }.set { ch_cellsam }

    //
    // Run the CELLSAM_SEGMENT_WBACKSUB subworkflow for samples that require
    // background subtraction and CellSAM segmentation
    //
    CELLSAM_SEGMENT_WBACKSUB(
        ch_cellsam.backsub_cellsam
    )

    //
    // Run CELLSAM_SEGMENT subworkflow for samples that ONLY require CellSAM segmentation
    //
    CELLSAM_SEGMENT(
        ch_cellsam.cellsam_only
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  + 'pipeline_software_' +  ''  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_collated_versions     // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
