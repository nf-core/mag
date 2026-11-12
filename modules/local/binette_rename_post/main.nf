process RENAME_POSTBINETTE {
    tag "${meta.assembler}-${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0'
        : 'biocontainers/seqkit:2.9.0--h9ee0642_0'}"

    input:
    tuple val(meta), path(bins), path(contig2bin), path(contigs)

    output:
    tuple val(meta), path("${meta.assembler}-BinetteRefined-${meta.id}.*.fa.gz", arity: '1..*'), emit: refined_bins
    tuple val(meta), path("${meta.assembler}-BinetteUnbinned-${meta.id}.fa.gz"), optional: true, emit: refined_unbins
    path "versions.yml", emit: versions

    script:
    def unbinned = "${meta.assembler}-BinetteUnbinned-${meta.id}.fa.gz"
    """
    # Binette names bins <prefix>_bin<N>.fa.gz, the rest of the pipeline expects <prefix>.<N>.fa.gz
    for bin in ${bins}; do
        num=\${bin##*_bin}
        num=\${num%.fa.gz}
        mv \${bin} ${meta.assembler}-BinetteRefined-${meta.id}.\${num}.fa.gz
    done

    # Binette does not report unbinned contigs, so take the contigs absent from every refined bin
    cut -f1 ${contig2bin} > binned_contigs.txt
    seqkit grep -v -f binned_contigs.txt ${contigs} -o ${unbinned}
    if [ -z "\$(seqkit seq -n ${unbinned} | head -n 1)" ]; then
        rm ${unbinned}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS
    """

    stub:
    """
    echo "" | gzip > ${meta.assembler}-BinetteRefined-${meta.id}.1.fa.gz
    echo "" | gzip > ${meta.assembler}-BinetteUnbinned-${meta.id}.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS
    """
}
