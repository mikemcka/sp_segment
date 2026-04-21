process SMOOTHMASKS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python_tifffile_numpy_scipy_pruned:e54488103afb8110"

    input:
    tuple val(meta), path(mask)

    output:
    tuple val(meta), path("*_smoothed.tiff"), emit: smoothed_mask
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_name = mask.getName()
    def output_name = input_name.replaceAll(/\.(tiff?|TIFF?)$/, '_smoothed.tiff')
    """
    smooth_masks.py \\
        ${mask} \\
        ${output_name} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        scipy: \$(python3 -c "import scipy; print(scipy.__version__)")
        scikit-image: \$(python3 -c "import skimage; print(skimage.__version__)")
        shapely: \$(python3 -c "import shapely; print(shapely.__version__)")
    END_VERSIONS
    """

    stub:
    def input_name = mask.getName()
    def output_name = input_name.replaceAll(/\.(tiff?|TIFF?)$/, '_smoothed.tiff')
    """
    touch ${output_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        scipy: \$(python3 -c "import scipy; print(scipy.__version__)")
        scikit-image: \$(python3 -c "import skimage; print(skimage.__version__)")
        shapely: \$(python3 -c "import shapely; print(shapely.__version__)")
    END_VERSIONS
    """
}
