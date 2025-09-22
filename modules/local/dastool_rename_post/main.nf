process RENAME_POSTDASTOOL {
    tag "${meta.assembler}-${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c262fc09eca59edb5a724080eeceb00fb06396f510aefb229c2d2c6897e63975/data'
        : 'community.wave.seqera.io/library/coreutils:9.5--ae99c88a9b28c264'}"

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path("${meta.assembler}-*Refined-${meta.id}.*.fa", includeInputs: true), optional: true, emit: refined_bins
    tuple val(meta), path("${meta.assembler}-DASToolUnbinned-${meta.id}.fa"), optional: true, emit: refined_unbins
    path "versions.yml", emit: versions

    script:
    def VERSION = '9.4.3'
    """
    if [[ -f unbinned.fa ]]; then
        mv unbinned.fa ${meta.assembler}-DASToolUnbinned-${meta.id}.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: ${VERSION}
    END_VERSIONS
    """
}
