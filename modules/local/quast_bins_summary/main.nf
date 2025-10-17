process QUAST_BINS_SUMMARY {

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    path summaries

    output:
    path ("quast_summary.tsv"), emit: summary
    path "versions.yml", emit: versions

    script:
    """
    QUAST_BIN=\$(echo \"${summaries}\" | sed 's/[][]//g')
    IFS=', ' read -r -a quast_bin <<< \"\$QUAST_BIN\"
    for quast_file in \"\${quast_bin[@]}\"; do
        if ! [ -f "quast_summary.tsv" ]; then
            cp "\${quast_file}" "quast_summary.tsv"
        else
            tail -n +2 "\${quast_file}" >> "quast_summary.tsv"
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cp: \$(cp --version 2>&1 | sed -n 1p | sed 's/cp (GNU coreutils) //')
    END_VERSIONS
    """
}
