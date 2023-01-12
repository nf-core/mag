process RENAME_PREDASTOOL {
    tag "${meta.assembler}-${meta.binner}-${meta.id}"
    label 'process_low'

    // Using container from multiqc since it'll be included anyway
    conda "bioconda::multiqc=1.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path("${meta.assembler}-${meta.binner}Refined-${meta.id}*"), emit: renamed_bins

    script:
    """
    if [ -n "${bins}" ]
    then
        for bin in ${bins}; do
            if [[ \${bin} =~ ${meta.assembler}-${meta.binner}-${meta.id}.([_[:alnum:]]+).fa ]]; then
                num=\${BASH_REMATCH[1]}
                mv \${bin} ${meta.assembler}-${meta.binner}Refined-${meta.id}.\${num}.fa
            else
                echo "ERROR: the bin filename \${bin} does not match the expected format '${meta.assembler}-${meta.binner}-${meta.id}.([_[:alnum:]]+).fa'!"
                exit 1
            fi
        done
    fi
    """
}
