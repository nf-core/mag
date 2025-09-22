process RENAME_PREDASTOOL {
    tag "${meta.assembler}-${meta.binner}-${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c262fc09eca59edb5a724080eeceb00fb06396f510aefb229c2d2c6897e63975/data'
        : 'community.wave.seqera.io/library/coreutils:9.5--ae99c88a9b28c264'}"

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path("${meta.assembler}-${meta.binner}Refined-${meta.id}*"), emit: renamed_bins
    path "versions.yml", emit: versions

    script:
    def VERSION = '9.4.3'
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: ${VERSION}
    END_VERSIONS
    """
}
