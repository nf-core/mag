process RENAME_POSTBINETTE {
    tag "${meta.assembler}-${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path("${meta.assembler}-BinetteRefined-${meta.id}.*.fa.gz", arity: '1..*'), emit: refined_bins
    path "versions.yml", emit: versions

    script:
    """
    # Binette names bins <prefix>_bin<N>.fa.gz, the rest of the pipeline expects <prefix>.<N>.fa.gz
    for bin in ${bins}; do
        num=\${bin##*_bin}
        num=\${num%.fa.gz}
        mv \${bin} ${meta.assembler}-BinetteRefined-${meta.id}.\${num}.fa.gz
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(echo \$(mv --version 2>&1) | sed 's/^.*(GNU coreutils) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    """
    echo "" | gzip > ${meta.assembler}-BinetteRefined-${meta.id}.1.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(echo \$(mv --version 2>&1) | sed 's/^.*(GNU coreutils) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
