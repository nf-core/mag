process DEEPMASED {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deepmased:0.3.1--pyh5ca1d4c_0':
        'biocontainers/deepmased:0.3.1--pyh5ca1d4c_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta)

    output:
    tuple val(meta), path("*_deepmased_predictions.tsv"), emit: predictions
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.3.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    echo -e "bam\\tfasta" > ${prefix}_file_paths.tsv
    echo -e "${bam}\\t${fasta}" >> ${prefix}_file_paths.tsv

    # Run features with -p 1 to avoid Python multiprocessing issues inside Docker
    DeepMAsED features \\
        ${prefix}_file_paths.tsv \\
        -p 1 \\
        -o . \\
        -n ${prefix}_feature_file_paths.tsv

    DeepMAsED predict \\
        ${prefix}_feature_file_paths.tsv \\
        --n-procs ${task.cpus} \\
        --cpu-only \\
        --save-name ${prefix}_deepmased \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepmased: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.3.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}_deepmased_predictions.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepmased: $VERSION
    END_VERSIONS
    """
}
