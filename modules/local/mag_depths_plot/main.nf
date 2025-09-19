process MAG_DEPTHS_PLOT {
    tag "${meta.assembler}-${meta.binner}-${meta.id}"
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-d14219255233ee6cacc427e28a7caf8ee42e8c91:0a22c7568e4a509925048454dad9ab37fa8fe776-0'
        : 'biocontainers/mulled-v2-d14219255233ee6cacc427e28a7caf8ee42e8c91:0a22c7568e4a509925048454dad9ab37fa8fe776-0'}"

    input:
    tuple val(meta), path(depths)
    path sample_groups

    output:
    tuple val(meta), path("${meta.assembler}-${meta.binner}-${meta.id}-binDepths.heatmap.png"), emit: heatmap
    path "versions.yml", emit: versions

    script:
    """
    plot_mag_depths.py --bin_depths ${depths} \
                    --groups ${sample_groups} \
                    --out "${meta.assembler}-${meta.binner}-${meta.id}-binDepths.heatmap.png"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
        seaborn: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('seaborn').version)")
    END_VERSIONS
    """
}
