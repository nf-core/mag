process METABINNER_BINS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabinner:1.4.4--hdfd78af_0' :
        'quay.io/biocontainers/metabinner:1.4.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(membership)

    output:
    tuple val(meta), path("*.tooShort.fa.gz")         , emit: tooshort
    tuple val(meta), path("*.unbinned.fa.gz")         , emit: unbinned
    tuple val(meta), path("bins/*.fa.gz")             , emit: bins
    path "versions.yml"                               , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_contig_size = task.ext.min_contig_size ?: "1000"
    """
    # unzip membership file
    zcat ${membership} > membership.tsv

    # collect bins & un-binned fractions
    create_metabinner_bins.py \\
        membership.tsv \\
        ${fasta} \\
        ./bins \\
        ${prefix} \\
        ${min_contig_size}
    find ./bins/ -name "*.fa" -type f | xargs -t -n 1 bgzip -@ ${task.cpus}

    # zip contig fractions
    find ./bins/ -name "*[tooShort,unbinned].fa.gz" -type f -exec mv {} . \\;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
