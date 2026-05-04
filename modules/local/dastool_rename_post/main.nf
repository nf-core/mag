process RENAME_POSTDASTOOL {
    tag "${meta.assembler}-${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

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
        coreutils: \$(echo \$(mv --version 2>&1) | sed 's/^.*(GNU coreutils) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
