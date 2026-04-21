/*
 * This subworkflow uses code adapted from nf-core sopa
 * Original source: https://github.com/nf-core/sopa
 * License: MIT
 */
include { SOPA_CELLPOSE      } from '../sopa_cellpose/main.nf'
include { SOPA_CONVERT       } from '../../../modules/local/sopa/convert/main.nf'
include { SOPA_PATCHIFYIMAGE } from '../../../modules/local/sopa/patchifyimage/main.nf'
include { PARQUETTOTIFF      } from '../../../modules/local/parquettotiff/main.nf'

workflow SOPA_SEGMENT_COMPARTMENT {

    take:
    ch_sopa // channel: [ (meta, tiff, nuclear_channel, membrane_channels) ]
    compartment

    main:

    ch_versions = channel.empty()

    ch_sopa.map {
        meta,
        tiff,
        _nuclear_channel,
        _membrane_channels -> [ meta, tiff ]
    }.set { ch_convert }

    //
    // Run SOPA convert to convert tiff to zarr format
    //
    SOPA_CONVERT(
        ch_convert,
        compartment
    )
    ch_versions = ch_versions.mix(SOPA_CONVERT.out.versions.first())


    //
    // Run SOPA patchify to create image patches
    //
    SOPA_PATCHIFYIMAGE(
        SOPA_CONVERT.out.spatial_data
    )
    ch_versions = ch_versions.mix(SOPA_PATCHIFYIMAGE.out.versions.first())

    // Create a channel for each patch
    SOPA_PATCHIFYIMAGE.out.patches
        .join( SOPA_CONVERT.out.spatial_data, by: 0 )
        .map { meta, patches_file_image, _image_patches, zarr ->
            [ meta, zarr, patches_file_image.text.trim().toInteger() ] }
        .flatMap { meta, zarr, n_patches ->
            (0..<n_patches).collect { index -> [ meta, zarr, index, n_patches ] } }
        .combine(ch_sopa, by: 0)
        .map { meta, zarr, index, n_patches, _tiff, nuclear_channel, membrane_channels ->
            [ meta, zarr, index, n_patches, nuclear_channel, membrane_channels ]
        }.set { ch_cellpose }

    //
    // Run SOPA with cellpose for nuclear segmentation
    //
    SOPA_CELLPOSE(
        ch_cellpose,
        SOPA_CONVERT.out.spatial_data
    )
    ch_versions = ch_versions.mix(SOPA_CELLPOSE.out.versions.first())

    //
    // Convert nuclear segmentation parquet to tiff
    //
    PARQUETTOTIFF(
        SOPA_CELLPOSE.out.boundaries
            .join(ch_sopa, by: 0)
            .map { meta, boundaries, tiff, _nuc_chan, _mem_chans ->
                [ meta, boundaries, tiff, compartment ]
            }
    )
    ch_versions = ch_versions.mix(PARQUETTOTIFF.out.versions.first())

    emit:
    zarr                 = SOPA_CONVERT.out.spatial_data         // channel: [ val(meta), *.zarr ]
    tiff                 = PARQUETTOTIFF.out.tiff                // channel: [ val(meta), *.tiff ]

    versions = ch_versions                                      // channel: [ versions.yml ]
}
