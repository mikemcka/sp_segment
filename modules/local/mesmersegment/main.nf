process MESMERSEGMENT {
    tag "$meta.id"
    label 'process_multi'
    secret 'DEEPCELL_ACCESS_TOKEN'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/wehi-soda-hub/mesmersegmentation:0.3.1'

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
    def membrane_channel_args = membrane_channels != '' && membrane_channels != [] ?
        membrane_channels.split(":").collect { channel ->
            "--membrane-channel \"${channel}\""
        }.join(' ') : ''
    """
    # Phase 1: Serialize model download — first task downloads, others wait
    REAL_HOME="${params.deepcell_cache_dir ? params.deepcell_cache_dir + '/$(whoami)' : '\$HOME'}"
    mkdir -p "\$REAL_HOME/.deepcell/models"
    (
        flock -x 200
        if ! ls "\$REAL_HOME/.deepcell/models"/MultiplexSegmentation*.tar.gz >/dev/null 2>&1; then
            HOME="\$REAL_HOME" python -c "from deepcell.applications import Mesmer; Mesmer()"
        fi
    ) 200>>"\$REAL_HOME/.deepcell/model_download.lock"

    # Phase 2: Task-local HOME so each task extracts the model independently
    # (deepcell always re-extracts the archive, causing corruption if concurrent)
    export HOME="\$PWD/.task_home"
    mkdir -p "\$HOME/.deepcell/models"
    for f in "\$REAL_HOME/.deepcell/models"/MultiplexSegmentation*.tar.gz; do
        [ -e "\$f" ] && ln -sf "\$f" "\$HOME/.deepcell/models/"
    done

    # Phase 3: Run segmentation
    mesmer-segment \\
        "${tiff}" \\
        --compartment ${compartment} \\
        --nuclear-channel ${nuclear_channel} \\
        ${membrane_channel_args} \\
        ${args} \\
        > "${prefix}_${compartment}.tiff"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mesmersegmentation: v0.1.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.tiff"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mesmersegmentation: v0.1.0
    END_VERSIONS
    """
}
