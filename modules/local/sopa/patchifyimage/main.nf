/*
 * This module uses code adapted from nf-core sopa
 * Original source: https://github.com/nf-core/sopa
 * License: MIT
 */
process SOPA_PATCHIFYIMAGE {
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container
        ? 'docker://quentinblampey/sopa:2.1.11'
        : 'docker.io/quentinblampey/sopa:2.1.11'}"

    input:
    tuple val(meta), path(zarr)

    output:
    tuple val(meta), path("*.zarr/.sopa_cache/patches_file_image"), path("*.zarr/shapes/image_patches"), emit: patches
    path "versions.yml"                                                                                , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    sopa patchify image \\
        ${args} \\
        ${zarr} \\
        --patch-width-pixel ${params.patch_width_pixel} \\
        --patch-overlap-pixel ${params.patch_overlap_pixel}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sopa: \$(sopa --version | sed 's/sopa //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.zarr/.sopa_cache
    mkdir -p ${prefix}.zarr/shapes

    touch ${prefix}.zarr/.sopa_cache/patches_file_image
    touch ${prefix}.zarr/shapes/image_patches

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sopa: \$(sopa --version | sed 's/sopa //')
    END_VERSIONS
    """
}
