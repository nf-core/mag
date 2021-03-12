// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process POOL_PAIRED_READS {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("pooled_${meta.id}_*.fastq.gz")

    script:
    """
    cat ${reads1} > "pooled_${meta.id}_R1.fastq.gz"
    cat ${reads2} > "pooled_${meta.id}_R2.fastq.gz"
    """
}
