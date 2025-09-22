process ADJUST_MAXBIN2_EXT {
    tag "${meta.assembler}-${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c262fc09eca59edb5a724080eeceb00fb06396f510aefb229c2d2c6897e63975/data'
        : 'community.wave.seqera.io/library/coreutils:9.5--ae99c88a9b28c264'}"

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path("*.fa.gz"), emit: renamed_bins
    path "versions.yml", emit: versions

    script:
    def VERSION = '9.4.3'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    if [ -n "${bins}" ]
    then
        for file in ${bins}; do
            [[ \${file} =~ (.*).fasta.gz ]];
            bin="\${BASH_REMATCH[1]}"
            mv \${file} \${bin}.fa.gz
        done
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: ${VERSION}
    END_VERSIONS
    """
}
