// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

/*
 * Raname FASTQ files using sample name to retrieve unique basenames
 */
process RENAME_FASTQS {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}{_R1,_R2,}.fastq.gz", includeInputs: true)

    script:
    if (!meta.single_end)
        """
        if ! [ -f "${meta.id}_R1.fastq.gz" ]; then
            ln -s "${reads[0]}" "${meta.id}_R1.fastq.gz"
        fi
        if ! [ -f "${meta.id}_R2.fastq.gz" ]; then
            ln -s "${reads[1]}" "${meta.id}_R2.fastq.gz"
        fi
        """
    else
        """
        if ! [ -f "${meta.id}.fastq.gz" ]; then
            ln -s "${reads}" "${meta.id}.fastq.gz"
        fi
        """
}
