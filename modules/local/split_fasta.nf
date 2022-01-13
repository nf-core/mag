process SPLIT_FASTA {
    tag "${meta.assembler}-${meta.id}"
    label 'process_low'

    // Using container from metabat2 process, since this will be anyway already downloaded and contains biopython and pandas
    conda (params.enable_conda ? "bioconda::metabat2=2.15 conda-forge::python=3.6.7 conda-forge::biopython=1.74 conda-forge::pandas=1.1.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-e25d1fa2bb6cbacd47a4f8b2308bd01ba38c5dd7:75310f02364a762e6ba5206fcd11d7529534ed6e-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-e25d1fa2bb6cbacd47a4f8b2308bd01ba38c5dd7:75310f02364a762e6ba5206fcd11d7529534ed6e-0"
    }

    input:
    tuple val(meta), path(unbinned)

    output:
    path "unbinned/*" , emit: unbinned

    script:
    """
    # save unbinned contigs above thresholds into individual files, dump others in one file
    split_fasta.py $unbinned ${params.min_length_unbinned_contigs} ${params.max_unbinned_contigs} ${params.min_contig_size}

    mkdir -p unbinned/
    mv *.fa unbinned/
    gzip unbinned/*
    """
}
