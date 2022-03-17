process SPLIT_FASTA {
    tag "${meta.assembler}-${meta.binner}-${meta.id}"
    label 'process_low'

    // Using container from metabat2 process, since this will be anyway already downloaded and contains biopython and pandas
    conda (params.enable_conda ? "bioconda::metabat2=2.15 conda-forge::python=3.6.7 conda-forge::biopython=1.74 conda-forge::pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e25d1fa2bb6cbacd47a4f8b2308bd01ba38c5dd7:75310f02364a762e6ba5206fcd11d7529534ed6e-0' :
        'quay.io/biocontainers/mulled-v2-e25d1fa2bb6cbacd47a4f8b2308bd01ba38c5dd7:75310f02364a762e6ba5206fcd11d7529534ed6e-0' }"

    input:
    tuple val(meta), path(unbinned)

    output:
    tuple val(meta), path("*.[0-9]*.fa.gz")     , optional:true, emit: unbinned
    tuple val(meta), path("*.pooled.fa.gz")     , optional:true, emit: pooled
    tuple val(meta), path("*.remaining.fa.gz")  , optional:true, emit: remaining
    path "versions.yml"                                 , emit: versions

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
