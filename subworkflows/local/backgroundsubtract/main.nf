include { EXTRACTMARKERS      } from '../../../modules/local/extractmarkers/main.nf'
include { BACKSUB             } from '../../../modules/nf-core/backsub/main.nf'

workflow BACKGROUNDSUBTRACT {

    take:
    ch_backsub // channel: background subtraction parameters (sample name and tiff)

    main:

    ch_versions = channel.empty()

    //
    // Extract markers from the input tiff file
    //
    EXTRACTMARKERS(
        ch_backsub
    )
    ch_versions = ch_versions.mix(EXTRACTMARKERS.out.versions.first())


    //
    // Run background subtraction module on tiff with extracted markers
    //
    BACKSUB(
        ch_backsub,
        EXTRACTMARKERS.out.markers
    )
    ch_versions = ch_versions.mix(BACKSUB.out.versions.first())

    emit:
    backsub_tif   = BACKSUB.out.backsub_tif    // channel: [ val(meta), *.ome.tif ]
    markers       = BACKSUB.out.markerout      // channel: [ val(meta2), markers.csv ]

    versions      = ch_versions                     // channel: [ versions.yml ]
}
