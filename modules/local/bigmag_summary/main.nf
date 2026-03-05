process PREPARE_BIGMAG_SUMMARY {

    conda "conda-forge::pandas=1.4.3"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pandas:1.4.3'
        : 'biocontainers/pandas:1.4.3'}"

    input:
    path summary
    path gunc_sum

    output:
    path "bigmag_summary.tsv", emit: bigmag_summary
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    summary  = summary.sort().size() > 0 ? "--summary ${summary}" : ""
    def gunc_summary  = gunc_sum.sort().size() > 0 ? "--gunc_summary ${gunc_sum}" : ""
    """
    prepare_bigmag_summary.py \
        ${args} \
        ${summary} \
        ${gunc_summary} \
        --out bigmag_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
