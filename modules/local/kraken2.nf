// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process KRAKEN2 {
    tag "${meta.id}-${db_name}"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::kraken2=2.0.8_beta" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/kraken2:2.0.8_beta--pl526hc9558a2_2"
    } else {
        container "quay.io/biocontainers/kraken2:2.0.8_beta--pl526hc9558a2_2"
    }

    input:
    tuple val(meta), path(reads)
    tuple val(db_name), path("database/*")

    output:
    tuple val("kraken2"), val(meta), path("results.krona"), emit: results_for_krona
    path  "kraken2_report.txt"
    path  '*.version.txt'                                  , emit: version

    script:
    def software = getSoftwareName(task.process)
    def input = meta.single_end ? "\"${reads}\"" :  "--paired \"${reads[0]}\" \"${reads[1]}\""
    """
    kraken2 \
        --report-zero-counts \
        --threads ${task.cpus} \
        --db database \
        --report kraken2_report.txt \
        $input \
        > kraken2.kraken
    cat kraken2.kraken | cut -f 2,3 > results.krona

    echo \$(kraken2 --version 2>&1) | sed 's/Kraken version //; s/ Copyright.*//' > ${software}.version.txt
    """
}