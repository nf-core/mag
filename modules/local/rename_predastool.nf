process RENAME_PREDASTOOL {
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
    tuple val(meta), path("${meta.assembler}-${meta.binner}Refined-${meta.id}*")           , optional:true, emit: renamed_bins

    script:
    """
    for i in ${bins}; do
        filename=\$i
        suffix=\${filename#*.}

        mv \$i ${meta.assembler}-${meta.binner}Refined-${meta.id}.\$suffix
    done
    """
}
