process METABINNER_METABINNER {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metabinner:1.4.4--hdfd78af_0' :
        'quay.io/biocontainers/metabinner:1.4.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(kmer), path(depth)
    val val_min_contig_size

    output:
    tuple val(meta), path("*.tsv.gz") , emit: membership
    tuple val(meta), path("*.log.gz") , emit: log
    path "versions.yml"               , emit: versions

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_contig_size = val_min_contig_size ?: "1000"
    def VERSION = '1.4.4-0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    # unzip kmer file
    zcat $kmer > kmer_profile.csv

    # create coverage profile in Metabinner format
    zcat ${depth} | awk '{if (\$2>${min_contig_size}) print \$0 }' | cut -f -1,4- > coverage_profile.tsv

    # requires absolute paths
    wd=\$(pwd)
    metabinner_path=\$(dirname \$(which run_metabinner.sh))
    bash run_metabinner.sh \\
        -a \${wd}/${fasta} \\
        -o \${wd}/${prefix} \\
        -d \${wd}/coverage_profile.tsv \\
        -k \${wd}/kmer_profile.csv \\
        -t ${task.cpus} \\
        -p \${metabinner_path} \\
        ${args}

    # clean-up unzipped kmer file
    rm kmer_profile.csv

    # collect & zip membership & log files
    gzip -cn ${prefix}/metabinner_res/metabinner_result.tsv > ${prefix}.tsv.gz
    gzip -cn ${prefix}/metabinner_res/result.log > ${prefix}.metabinner.log.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MetaBinner: $VERSION
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
