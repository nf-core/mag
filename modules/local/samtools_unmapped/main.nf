process SAMTOOLS_UNMAPPED {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9e/9edc2564215d5cd137a8b25ca8a311600987186d406b092022444adf3c4447f7/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.21--6cb89bfd40cbaabf'}"

    input:
    tuple val(meta), path(input), path(index)

    output:
    tuple val(meta), path("*hostremoved.fastq.gz"), emit: fastq
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mapped = "-o ${prefix}_mapped.fastq.gz"

    """
    samtools \\
        view \\
        --threads ${task.cpus - 1} \\
        ${args} \\
        ${input} \\
        | \\
    samtools \\
        fastq \\
        ${args2} \\
        --threads ${task.cpus - 1} \\
        -0 ${prefix}.fastq.gz \\
        ${mapped}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
