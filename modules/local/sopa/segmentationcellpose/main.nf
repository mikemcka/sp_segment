/*
 * This module uses code adapted from nf-core sopa
 * Original source: https://github.com/nf-core/sopa
 * License: MIT
 */
process SOPA_SEGMENTATIONCELLPOSE {
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container
        ? 'docker://quentinblampey/sopa:2.1.11-cellpose'
        : 'docker.io/quentinblampey/sopa:2.1.11-cellpose'}"

    input:
    tuple val(meta), path(zarr), val(index), val(n_patches), val(nuclear_channel), val(membrane_channel)

    output:
    tuple val(meta), path("*.zarr/.sopa_cache/cellpose_boundaries/${index}.parquet"), emit: cellpose_parquet
    path "versions.yml"                                                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    def membrane_channel_arg = (membrane_channel && membrane_channel != "[]") ? "--channels \"${membrane_channel}\"" : ""
    """
    export NUMBA_CACHE_DIR=\$PWD/.numba_cache

    sopa segmentation cellpose \\
        ${args} \\
        --patch-index ${index} \\
        ${membrane_channel_arg} \\
        --channels "${nuclear_channel}" \\
        --diameter ${params.cellpose_diameter} \\
        --min-area ${params.cellpose_min_area} \\
        ${zarr}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sopa: \$(sopa --version | sed 's/sopa //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.zarr/.sopa_cache/cellpose_boundaries
    touch ${prefix}.zarr/.sopa_cache/cellpose_boundaries/${index}.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sopa: \$(sopa --version | sed 's/sopa //')
    END_VERSIONS
    """
}
