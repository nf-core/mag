process MAG_DEPTHS_SUMMARY {

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.4.3' :
        'biocontainers/pandas:1.4.3' }"

    input:
    path(mag_depths)

    output:
    path("${prefix}.tsv"), emit: summary
    path "versions.yml"  , emit: versions

    script:
    prefix = task.ext.prefix ?: "bin_depths_summary"
    """
    get_mag_depths_summary.py --depths ${mag_depths} \
                            --out "${prefix}.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
