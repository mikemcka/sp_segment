process CELLMEASUREMENT {
    tag "$meta.id"
    label 'process_multi'

    conda "${moduleDir}/environment.yml"
    container 'ghcr.io/wehi-soda-hub/cellmeasurement-py:0.1.3'

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
    def nuclear_mask_flag = params.use_whole_cell_only ? '' : "--nuclear-mask ${nuclear_mask}"
    """
    cellmeasurement \\
        --whole-cell-mask ${whole_cell_mask} \\
        ${nuclear_mask_flag} \\
        --tiff-file ${tiff} \\
        --output-file ${prefix}.geojson \\
        --output-mask ${prefix}_mask.tiff \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellmeasurement: \$(cellmeasurement --version)
        python: \$(python3 --version | sed 's/Python //')
        scipy: \$(python3 -c "import scipy; print(scipy.__version__)")
        scikit-image: \$(python3 -c "import skimage; print(skimage.__version__)")
        shapely: \$(python3 -c "import shapely; print(shapely.__version__)")
        tifffile: \$(python3 -c "import tifffile; print(tifffile.__version__)")
        dask: \$(python3 -c "import dask; print(dask.__version__)")
        spatialdata: \$(python3 -c "import spatialdata; print(spatialdata.__version__)")
        rasterio: \$(python3 -c "import rasterio; print(rasterio.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        geopandas: \$(python3 -c "import geopandas; print(geopandas.__version__)")
        typer: \$(python3 -c "import typer; print(typer.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.geojson
    touch ${prefix}_mask.tiff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellmeasurement: \$(cellmeasurement --version)
        python: \$(python3 --version | sed 's/Python //')
        scipy: \$(python3 -c "import scipy; print(scipy.__version__)")
        scikit-image: \$(python3 -c "import skimage; print(skimage.__version__)")
        shapely: \$(python3 -c "import shapely; print(shapely.__version__)")
        tifffile: \$(python3 -c "import tifffile; print(tifffile.__version__)")
        dask: \$(python3 -c "import dask; print(dask.__version__)")
        spatialdata: \$(python3 -c "import spatialdata; print(spatialdata.__version__)")
        rasterio: \$(python3 -c "import rasterio; print(rasterio.__version__)")
        numpy: \$(python3 -c "import numpy; print(numpy.__version__)")
        geopandas: \$(python3 -c "import geopandas; print(geopandas.__version__)")
        typer: \$(python3 -c "import typer; print(typer.__version__)")
    END_VERSIONS
    """
}
