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
    # use find as sometimes these are empty and need to fail gracefully
    find -L -type f -name "*bin2classification.names.txt.gz" -exec sh -c 'for f do gunzip -c \$f > \${f%.*}; done' find-sh {} +

    bioawk '(NR == 1) || (FNR > 1)' *.txt > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioawk: \$(bioawk --version | cut -f 3 -d ' ' )
    END_VERSIONS
    """
}
