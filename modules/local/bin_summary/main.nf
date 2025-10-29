process BIN_SUMMARY {

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:1.4.3'
        : 'biocontainers/pandas:1.4.3'}"

    input:
    path bin_depths
    path quast_sum
    path gtdbtk_sum
    path cat_sum
    path busco_sum
    path checkm_sum
    path checkm2_sum

    output:
    path "bin_summary.tsv", emit: summary
    path "versions.yml", emit: versions

    script:
    def quast_summary = quast_sum.sort().size() > 0 ? "--quast_summary ${quast_sum}" : ""
    def gtdbtk_summary = gtdbtk_sum.sort().size() > 0 ? "--gtdbtk_summary ${gtdbtk_sum}" : ""
    def cat_summary = cat_sum.sort().size() > 0 ? "--cat_summary ${cat_sum}" : ""
    def busco_summary = busco_sum.sort().size() > 0 ? "--busco_summary ${busco_sum}" : ""
    def checkm_summary = checkm_sum.sort().size() > 0 ? "--checkm_summary ${checkm_sum}" : ""
    def checkm2_summary = checkm2_sum.sort().size() > 0 ? "--checkm2_summary ${checkm2_sum}" : ""
    """
    combine_tables.py \
        --depths_summary ${bin_depths} \
        ${quast_summary} \
        ${gtdbtk_summary} \
        ${cat_summary} \
        ${busco_summary} \
        ${checkm_summary} \
        ${checkm2_summary} \
        --out bin_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
