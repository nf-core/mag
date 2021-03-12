// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process POOL_SINGLE_READS {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("pooled_${meta.id}.fastq.gz")

    script:
    """
    cat ${reads} > "pooled_${meta.id}.fastq.gz"
    """
}
