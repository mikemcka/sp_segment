process SEGMENTATIONREPORT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/wehi-soda-hub/spatialvis:0.2.0'

    input:
    tuple val(meta),
        path(annotations),
        val(run_mesmer),
        val(run_cellpose),
        val(run_cellsam),
        val(nuclear_channel),
        val(membrane_channels),
        path(image_file)

    output:
    tuple val(meta), path("*/*.html")       , emit: report
    tuple val(meta), path("*/*.rds")        , emit: rds, optional: true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    # Process-specific cache directory to avoid cache conflicts
    export XDG_CACHE_HOME="\$(pwd)/.quarto_cache"
    mkdir -p "\$XDG_CACHE_HOME"

    Rscript -e "spatialVis::copy_report_template(
        template_name = 'segmentation_report_template.qmd',
        output_dir = '.',
        overwrite = TRUE
    )"
    quarto render segmentation_report_template.qmd \\
        --to html \\
        --no-cache \\
        --output ${prefix}.html \\
        ${args} \\
        -P geojson_file:${annotations} \\
        -P sample_name:${meta.id} \\
        -P nuclear_channel:${nuclear_channel} \\
        -P membrane_channels:"${membrane_channels}" \\
        -P image_file:"${image_file}" \\
        -P run_cellpose:${run_cellpose} \\
        -P run_mesmer:${run_mesmer} \\
        -P run_cellsam:${run_cellsam}

    mkdir -p ${prefix}
    mv ${prefix}.html ${prefix}
    mv segmentation_report_template.qmd ${prefix}/${prefix}.qmd
    if [[ -f ${prefix}.rds ]]; then
        mv ${prefix}.rds ${prefix}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
        spatialVis: \$(Rscript -e "cat(as.character(packageVersion('spatialVis')))")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}.html
    touch ${prefix}/${prefix}.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
        spatialVis: \$(Rscript -e "cat(as.character(packageVersion('spatialVis')))")
    END_VERSIONS
    """
}
