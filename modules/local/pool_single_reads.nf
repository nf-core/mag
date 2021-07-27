// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process POOL_SINGLE_READS {
    tag "$meta.id"

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("pooled_${meta.id}.fastq.gz"), emit: reads

    script:
    """
    cat ${reads} > "pooled_${meta.id}.fastq.gz"
    """
}
