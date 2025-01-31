process CAT_SUMMARY {
    label 'process_low'

    conda "bioconda::bioawk=1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioawk:1.0--hed695b0_5' :
        'biocontainers/bioawk:1.0--hed695b0_5' }"

    input:
    path(cat_summaries)

    output:
    path("*.tsv")      , emit: combined
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "cat_summary"
    """
    bioawk '(NR == 1) || (FNR > 1)' *.txt > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioawk: \$(bioawk --version | cut -f 3 -d ' ' )
    END_VERSIONS
    """
}
