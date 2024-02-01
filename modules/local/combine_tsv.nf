process COMBINE_TSV {

    // Using bioawk as already use that for CONVERT_DEPTHS and does same thing
    conda "bioconda::bioawk=1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioawk:1.0--hed695b0_5' :
        'biocontainers/bioawk:1.0--hed695b0_5' }"

    input:
    path(bin_summaries)

    output:
    path("*.tsv")      , emit: combined
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "bin_depths_summary_combined"
    """
    bioawk '(NR == 1) || (FNR > 1)' ${bin_summaries} > ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioawk: \$(bioawk --version | cut -f 3 -d ' ' )
    END_VERSIONS
    """
}
