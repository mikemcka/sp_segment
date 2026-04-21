process COMBINECHANNELS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/tifffile_xarray_numpy_typer:f92759840da2dc33'

    input:
    tuple val(meta), path(tiff), val(nuclear_channel), val(membrane_channels)

    output:
    tuple val(meta), path("${prefix}_combined_channels.tiff"), emit: combined_tiff
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def membrane_channel_args = membrane_channels != '' && membrane_channels != [] ?
        membrane_channels.split(":").collect { channel ->
            "--membrane-channel \"${channel}\""
        }.join(' ') : ''
    """
    combine_channels.py \\
        $args \\
        --nuclear-channel "${nuclear_channel}" \\
        ${membrane_channel_args} \\
        $tiff > ${prefix}_combined_channels.tiff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_combined_channels.tiff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
