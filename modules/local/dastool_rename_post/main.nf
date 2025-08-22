process RENAME_POSTDASTOOL {
    tag "${meta.assembler}-${meta.id}"
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:20.04'
        : 'nf-core/ubuntu:20.04'}"

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path("${meta.assembler}-*Refined-${meta.id}.*.fa", includeInputs: true), optional: true, emit: refined_bins
    tuple val(meta), path("${meta.assembler}-DASToolUnbinned-${meta.id}.fa"), optional: true, emit: refined_unbins
    path "versions.yml", emit: versions

    script:
    """
    if [[ -f unbinned.fa ]]; then
        mv unbinned.fa ${meta.assembler}-DASToolUnbinned-${meta.id}.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(apt-cache policy coreutils | sed '2!d; s/.* //')
    END_VERSIONS
    """
}
