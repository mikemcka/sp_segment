/*
 * This module uses code adapted from nf-core sopa
 * Original source: https://github.com/nf-core/sopa
 * License: MIT
 */
process SOPA_RESOLVECELLPOSE {
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container
        ? 'docker://quentinblampey/sopa:2.2.6'
        : 'docker.io/quentinblampey/sopa:2.2.6'}"

    input:
    tuple val(meta), path(zarr), val(cellpose_parquet)

    output:
    tuple val(meta), path("*.zarr/shapes/cellpose_boundaries/*.parquet"), emit: cellpose_boundaries
    path "versions.yml"                                                 , emit: versions

    script:
    """
    # The apptainer profile runs with --no-home, so numba (imported by sopa) has no
    # writable cache directory and fails with "no locator available".
    export NUMBA_CACHE_DIR=\$PWD/.numba_cache

    sopa resolve cellpose ${zarr}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sopa: \$(sopa --version | sed 's/sopa //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.zarr/shapes/cellpose_boundaries
    touch ${prefix}.zarr/shapes/cellpose_boundaries/resolved.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sopa: \$(sopa --version | sed 's/sopa //')
    END_VERSIONS
    """
}
