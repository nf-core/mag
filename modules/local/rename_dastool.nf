process RENAME_DASTOOL {
    tag "${meta.assembler}-${meta.id}"
    label 'process_low'

    // Using container from multiqc since it'll be included anyway
    conda (params.enable_conda ? "bioconda::multiqc=1.12" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path("*_refined.fasta")           , optional:true, emit: refined_bins
    tuple val(meta), path("*_refined_unbinned.fasta")  , optional:true, emit: refined_unbins

    script:
    """
    if [[ -f unbinned.fa ]]; then
        mv unbinned.fa ${meta.assembler}-${meta.id}_refined_unbinned.fasta
    fi

    for i in *.fa; do
        mv \$i \${i%%.fa}_refined.fasta
    done
    """
}
