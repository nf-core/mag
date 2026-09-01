process METABINNER_TOOSHORT {
    tag "$meta.id"
    label 'process_low'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabinner:1.4.4--hdfd78af_0' :
        'quay.io/biocontainers/metabinner:1.4.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    val val_min_contig_size

    output:
    tuple val(meta), path("*.${min_contig_size}.fa") , emit: sizefiltered
    path "versions.yml"                              , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    min_contig_size = val_min_contig_size ?: "1000"
    def VERSION = '1.4.4-0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    metabinner_path=\$(dirname \$(which run_metabinner.sh))

    # filter contigs < ${min_contig_size} bp from fasta (default 1000)
    python \${metabinner_path}/scripts/Filter_tooshort.py ${fasta} ${min_contig_size}

    mv ${fasta.baseName}_${min_contig_size}.fa ${prefix}.${min_contig_size}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MetaBinner: $VERSION
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
