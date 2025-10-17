process GTDBTK_SUMMARY {


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'biocontainers/pandas:1.4.3' }"

    input:
    path(qc_discarded_bins)
    path(gtdbtk_summaries)
    path(filtered_bins)
    path(failed_bins)

    output:
    path "gtdbtk_summary.tsv", emit: summary
    path "versions.yml"      , emit: versions

    script:
    def args = task.ext.args ?: ''
    def discarded = qc_discarded_bins.sort().size() > 0 ? "--qc_discarded_bins ${qc_discarded_bins}" : ""
    def summaries = gtdbtk_summaries.sort().size() > 0 ?  "--summaries ${gtdbtk_summaries}" : ""
    def filtered  = filtered_bins.sort().size() > 0 ?     "--filtered_bins ${filtered_bins}" : ""
    def failed    = failed_bins.sort().size() > 0 ?       "--failed_bins ${failed_bins}" : ""
    """
    summary_gtdbtk.py $args $discarded $summaries $filtered $failed --out gtdbtk_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
