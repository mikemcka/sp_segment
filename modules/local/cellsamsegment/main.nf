process CELLSAMSEGMENT {
    tag "$meta.id"
    label 'process_gpu'
    secret 'DEEPCELL_ACCESS_TOKEN'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/python_tifffile_scikit-image_scikit-learn_pruned:593e00ba324c12b3'

    input:
    tuple val(meta), path(tiff), val(nuclear_channel), val(membrane_channels)
    val(compartment)

    output:
    tuple val(meta), path("*.tiff"), emit: segmentation_mask
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_channels = membrane_channels != '' && membrane_channels != [] ? membrane_channels.split(":") : []
    def mem_channel_args = mem_channels.collect { channel ->
        "--membrane-channel '${channel}'"
        }.join(' ')
    def cache_home = params.deepcell_cache_dir ? "export HOME=\"${params.deepcell_cache_dir}/\$(whoami)\"" : ''

    """
    ${cache_home}
    cellsam_segment.py \\
        ${tiff} \\
        --output ${prefix}_${compartment}.tiff \\
        --compartment ${compartment} \\
        --nuclear-channel '${nuclear_channel}' \\
        ${mem_channel_args} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellsam: 0.1.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_${compartment}.tiff"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellsam: 0.1.0
    END_VERSIONS
    """
}
