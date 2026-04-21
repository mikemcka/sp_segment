process CELLMEASUREMENT {
    tag "$meta.id"
    label 'process_multi'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_tifffile_numpy_scipy_pruned:e54488103afb8110"

    input:
    tuple val(meta),
        path(tiff),
        path(nuclear_mask, stageAs: 'nuclear_mask_input.tiff'),
        path(whole_cell_mask, stageAs: 'whole_cell_mask_input.tiff')

    output:
    tuple val(meta), path("*.geojson{,.gz}"), emit: annotations
    tuple val(meta), path("*_mask.tiff"), emit: masks
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cellmeasurement.py \\
        --nuclear-mask ${nuclear_mask} \\
        --whole-cell-mask ${whole_cell_mask} \\
        --tiff-file ${tiff} \\
        --output-file ${prefix}.geojson \\
        --output-mask ${prefix}_mask.tiff \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        scipy: \$(python3 -c "import scipy; print(scipy.__version__)")
        scikit-image: \$(python3 -c "import skimage; print(skimage.__version__)")
        shapely: \$(python3 -c "import shapely; print(shapely.__version__)")
        tifffile: \$(python3 -c "import tifffile; print(tifffile.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.geojson
    touch ${prefix}_mask.tiff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        scipy: \$(python3 -c "import scipy; print(scipy.__version__)")
        scikit-image: \$(python3 -c "import skimage; print(skimage.__version__)")
        shapely: \$(python3 -c "import shapely; print(shapely.__version__)")
        tifffile: \$(python3 -c "import tifffile; print(tifffile.__version__)")
    END_VERSIONS
    """
}
