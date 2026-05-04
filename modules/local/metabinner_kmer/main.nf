process METABINNER_KMER {
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
    tuple val(meta), path("*_kmer_4_f${min_contig_size}.csv.gz"), emit: composition_profile
    path "versions.yml"                                         , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    min_contig_size = val_min_contig_size ?: "1000"
    def VERSION = '1.4.4-0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    metabinner_path=\$(dirname \$(which run_metabinner.sh))

    # create composition profile (contigs > ${min_contig_size} p (default 1000), k = 4)
    python \${metabinner_path}/scripts/gen_kmer.py ${fasta} ${min_contig_size} 4

    gzip -cn ${fasta.baseName}_kmer_4_f${min_contig_size}.csv > ${prefix}_kmer_4_f${min_contig_size}.csv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MetaBinner: $VERSION
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
