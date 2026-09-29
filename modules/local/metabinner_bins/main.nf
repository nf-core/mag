process METABINNER_BINS {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabinner:1.4.4--hdfd78af_0' :
        'quay.io/biocontainers/metabinner:1.4.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(membership)

    output:
    tuple val(meta), path("*.unbinned.fa.gz"),            emit: unbinned, optional: true
    tuple val(meta), path("bins/*.fa.gz", arity: '1..*'), emit: bins
    path "versions.yml",                                  emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # unzip membership file
    zcat ${membership} > membership.tsv

    create_metabinner_bins.py \\
        membership.tsv \\
        ${fasta} \\
        ./bins \\
        ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
