/*
 * This module uses code adapted from nf-core sopa
 * Original source: https://github.com/nf-core/sopa
 * License: MIT
 */
process SOPA_CONVERT {
    label "process_high"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container
        ? 'docker://quentinblampey/sopa:2.1.11'
        : 'docker.io/quentinblampey/sopa:2.1.11'}"

    input:
    tuple val(meta), path(tiff)
    val(compartment)

    output:
    tuple val(meta), path("*.zarr"), emit: spatial_data
    path "versions.yml"            , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    sopa convert \\
        ${args} \\
        --sdata-path ${meta.id}_${compartment}.zarr \\
        --technology ${params.technology} \\
        ${tiff}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sopa: \$(sopa --version | sed 's/sopa //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${compartment}.zarr

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sopa: \$(sopa --version | sed 's/sopa //')
    END_VERSIONS
    """
}
