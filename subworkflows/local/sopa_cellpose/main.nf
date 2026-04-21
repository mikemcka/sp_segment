/*
 * This subworkflow uses code adapted from nf-core sopa
 * Original source: https://github.com/nf-core/sopa
 * License: MIT
 */

include { SOPA_SEGMENTATIONCELLPOSE } from '../../../modules/local/sopa/segmentationcellpose/main.nf'
include { SOPA_RESOLVECELLPOSE      } from '../../../modules/local/sopa/resolvecellpose/main.nf'

workflow SOPA_CELLPOSE {

    take:
    ch_patches // channel: [ (meta, zarr, index, n_patches, nuclear_channel, membrane_channels) ]
    ch_spatial_data // channel: [ (meta, zarr) ]

    main:

    ch_versions = channel.empty()

    //
    // Run SOPA segmentation with cellpose
    //
    SOPA_SEGMENTATIONCELLPOSE(
        ch_patches
    )
    ch_versions = ch_versions.mix(SOPA_SEGMENTATIONCELLPOSE.out.versions.first())


    // Collect cellpose segmentation boundaries into one channel per sample
    SOPA_SEGMENTATIONCELLPOSE.out.cellpose_parquet
        .groupTuple()
        .join( ch_spatial_data, by: 0 )
        .map { meta, cellpose_parquet, zarr ->
            [ meta, zarr, cellpose_parquet ]
        }
        .set { ch_resolve_cellpose }

    //
    // Resolve Cellpose segmentation boundaries
    //
    SOPA_RESOLVECELLPOSE(
        ch_resolve_cellpose
    )
    ch_versions = ch_versions.mix(SOPA_RESOLVECELLPOSE.out.versions.first())

    emit:
    boundaries  = SOPA_RESOLVECELLPOSE.out.cellpose_boundaries  // channel: [ val(meta), *.zarr/shapes/cellose_boundaries/*.parquet ]

    versions = ch_versions                                     // channel: [ versions.yml ]
}
