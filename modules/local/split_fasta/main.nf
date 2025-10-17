process SPLIT_FASTA {
    tag "${meta.assembler}-${meta.binner}-${meta.id}"
    label 'process_low'

    // Using container from metabat2 process, since this will be anyway already downloaded and contains biopython and pandas
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e25d1fa2bb6cbacd47a4f8b2308bd01ba38c5dd7:75310f02364a762e6ba5206fcd11d7529534ed6e-0' :
        'biocontainers/mulled-v2-e25d1fa2bb6cbacd47a4f8b2308bd01ba38c5dd7:75310f02364a762e6ba5206fcd11d7529534ed6e-0' }"

    input:
    tuple val(meta), path(unbinned)

    output:
    tuple val(meta), path("${meta.assembler}-${meta.binner}-${meta.id}.*.[1-9]*.fa.gz")   , optional:true, emit: unbinned
    tuple val(meta), path("${meta.assembler}-${meta.binner}-${meta.id}.*.pooled.fa.gz")   , optional:true, emit: pooled
    tuple val(meta), path("${meta.assembler}-${meta.binner}-${meta.id}.*.remaining.fa.gz"), optional:true, emit: remaining
    path "versions.yml"                                                                   , emit: versions

    script:
    """
    # save unbinned contigs above thresholds into individual files, dump others in one file
    split_fasta.py $unbinned ${params.min_length_unbinned_contigs} ${params.max_unbinned_contigs} ${params.min_contig_size}

    gzip *.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        biopython: 1.7.4
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
