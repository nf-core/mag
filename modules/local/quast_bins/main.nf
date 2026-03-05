process QUAST_BINS {
    tag "${meta.assembler}-${meta.binner}-${meta.domain}-${meta.refinement}-${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.0.2--py37pl526hb5aa323_2' :
        'biocontainers/quast:5.0.2--py37pl526hb5aa323_2' }"

    input:
    tuple val(meta), path(bins)

    output:
    path "QUAST/*", type: 'dir'     , emit: dir
    tuple val(meta), path("QUAST/*-quast_summary.tsv"), emit: quast_bin_summaries
    path "versions.yml"             , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.assembler}-${meta.binner}-${meta.domain}-${meta.refinement}-${meta.id}"
    """
    BINS=\$(echo \"$bins\" | sed 's/[][]//g')
    IFS=', ' read -r -a bins <<< \"\$BINS\"
    for bin in \"\${bins[@]}\"; do
        metaquast.py --threads "${task.cpus}" --max-ref-number 0 --rna-finding --gene-finding -l "\${bin}" "\${bin}" -o "QUAST/\${bin}"
        if ! [ -f "QUAST/${prefix}-quast_summary.tsv" ]; then
            cp "QUAST/\${bin}/transposed_report.tsv" "QUAST/${prefix}-quast_summary.tsv"
        else
            tail -n +2 "QUAST/\${bin}/transposed_report.tsv" >> "QUAST/${prefix}-quast_summary.tsv"
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        metaquast: \$(metaquast.py --version | sed "s/QUAST v//; s/ (MetaQUAST mode)//")
    END_VERSIONS
    """
}
