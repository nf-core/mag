process ADJUST_MAXBIN2_EXT {
    tag "${meta.assembler}-${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path("*.fa.gz"), emit: renamed_bins
    path "versions.yml", emit: versions

    script:
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
        coreutils: \$(echo \$(mv --version 2>&1) | sed 's/^.*(GNU coreutils) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
