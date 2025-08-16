process BIN_SUMMARY {

    conda "conda-forge::pandas=1.4.3"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:1.4.3'
        : 'biocontainers/pandas:1.4.3'}"

    input:
    path bin_depths
    path binqc_sum
    path quast_sum
    path gtdbtk_sum
    path cat_sum
    val  binqc_tool

    output:
    path "bin_summary.tsv", emit: summary
    path "versions.yml"   , emit: versions

    script:
    def binqc_summary  = binqc_sum.sort().size() > 0 ? "--binqc_summary ${binqc_sum}" : ""
    def quast_summary  = quast_sum.sort().size() > 0 ? "--quast_summary ${quast_sum}" : ""
    def gtdbtk_summary = gtdbtk_sum.sort().size() > 0 ? "--gtdbtk_summary ${gtdbtk_sum}" : ""
    def cat_summary    = cat_sum.sort().size() > 0 ? "--cat_summary ${cat_sum}" : ""
    """
    combine_tables.py \
        --depths_summary ${bin_depths} \
        ${binqc_summary} \
        ${quast_summary} \
        ${gtdbtk_summary} \
        ${cat_summary} \
        --binqc_tool ${binqc_tool} \
        --out bin_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
