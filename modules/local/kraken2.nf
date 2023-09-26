process KRAKEN2 {
    tag "${meta.id}-${db_name}"

    conda "bioconda::kraken2=2.0.8_beta"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kraken2:2.0.8_beta--pl526hc9558a2_2' :
        'biocontainers/kraken2:2.0.8_beta--pl526hc9558a2_2' }"

    input:
    tuple val(meta), path(reads)
    tuple val(db_name), path("database/*")

    output:
    tuple val("kraken2"), val(meta), path("results.krona"), emit: results_for_krona
    tuple val(meta), path("*kraken2_report.txt")          , emit: report
    path "versions.yml"                                   , emit: versions

    script:
    def input = meta.single_end ? "\"${reads}\"" :  "--paired \"${reads[0]}\" \"${reads[1]}\""
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    kraken2 \
        --report-zero-counts \
        --threads ${task.cpus} \
        --db database \
        --report ${prefix}.kraken2_report.txt \
        $input \
        > kraken2.kraken
    cat kraken2.kraken | cut -f 2,3 > results.krona

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //' | sed 's/ Copyright.*//')
    END_VERSIONS
    """
}
