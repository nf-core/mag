process SUMMARISE_PYDAMAGEBINS {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:1.4.3'
        : 'biocontainers/pandas:1.4.3'}"

    input:
    path pydamage_reports
    path contig_to_bin_map

    output:
    path "pydamage_bins_summary.tsv", emit: pydamage_bin_summary
    path "*_pydamage_bin_results.tsv", emit: pydamage_bin_results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    summarise_pydamagebins.py \\
        ${args} \\
        --contig-to-bin-map ${contig_to_bin_map} \\
        --output pydamage_bins_summary.tsv \\
        --verbose \\
        ${pydamage_reports.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    echo ${args}
    touch pydamage_bins_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
