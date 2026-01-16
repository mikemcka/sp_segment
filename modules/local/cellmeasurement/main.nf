process CELLMEASUREMENT {
    tag "$meta.id"
    label 'process_multi'

    conda "${moduleDir}/environment.yml"
    container "ghcr.io/wehi-soda-hub/cellmeasurement:0.2.3"

    input:
    tuple val(meta),
        path(tiff),
        path(nuclear_mask),
        path(whole_cell_mask)

    output:
    tuple val(meta), path("*.geojson"), emit: annotations
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    /cellmeasurement.sh \\
        --args="${args} \\
            --nuclear-mask=\$(readlink ${nuclear_mask}) \\
            --whole-cell-mask=\$(readlink ${whole_cell_mask}) \\
            --tiff-file=\$(readlink ${tiff}) \\
            --output-file=\$PWD/${prefix}.geojson \\
            --threads=${task.cpus}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellmeasurement: \$(/cellmeasurement.sh --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.geojson

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellmeasurement: \$(/cellmeasurement.sh --version)
    END_VERSIONS
    """
}
