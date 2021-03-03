// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process FILTLONG {
    tag "$meta.id"

    conda (params.enable_conda ? "bioconda::filtlong=0.2.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/filtlong:0.2.0--he513fc3_3"
    } else {
        container "quay.io/biocontainers/filtlong:0.2.0--he513fc3_3"
    }

    input:
    tuple val(meta), path(long_reads), path(short_reads_1), path(short_reads_2)

    output:
    tuple val(meta), path("${meta.id}_lr_filtlong.fastq.gz"), emit: reads
    path '*.version.txt'                                    , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    filtlong \
        -1 ${short_reads_1} \
        -2 ${short_reads_2} \
        --min_length ${params.longreads_min_length} \
        --keep_percent ${params.longreads_keep_percent} \
        --trim \
        --length_weight ${params.longreads_length_weight} \
        ${long_reads} | gzip > ${meta.id}_lr_filtlong.fastq.gz

    filtlong --version | sed -e "s/Filtlong v//g" > ${software}.version.txt
    """
}

