process BUSCO_SUMMARY {

    conda "conda-forge::pandas=1.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'biocontainers/pandas:1.4.3' }"

    input:
    path(summaries_domain)
    path(summaries_specific)
    path(failed_bins)

    output:
    path "busco_summary.tsv", emit: summary
    path "versions.yml"     , emit: versions

    script:
    def reference = params.busco_db.toString().contains('odb10')
    def auto = reference ? "" : "-a"
    def ss = summaries_specific.sort().size() > 0 ? "-ss ${summaries_specific}" : ""
    def sd = summaries_domain.sort().size() > 0 ? "-sd ${summaries_domain}" : ""
    def f = ""
    if ("${reference}" == false && failed_bins.sort().size() > 0)
        f = "-f ${failed_bins}"
    """
    summary_busco.py $auto $ss $sd $f -o busco_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}

