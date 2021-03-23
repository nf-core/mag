// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process POOL_PAIRED_READS {
    tag "$meta.id"

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("pooled_${meta.id}_*.fastq.gz"), emit: reads

    script:
    """
    cat ${reads1} > "pooled_${meta.id}_1.fastq.gz"
    cat ${reads2} > "pooled_${meta.id}_2.fastq.gz"
    """
}
