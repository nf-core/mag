process POOL_PAIRED_READS {
    label 'process_single'
    tag "${meta.id}"
    conda "conda-forge::sed=4.7"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:20.04'
        : 'nf-core/ubuntu:20.04'}"

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("pooled_${meta.id}_*.fastq.gz"), emit: reads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat ${reads1} > "pooled_${meta.id}_1.fastq.gz"
    cat ${reads2} > "pooled_${meta.id}_2.fastq.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(cat --version 2>&1 | sed -n 1p | sed 's/cat (GNU coreutils) //')
    END_VERSIONS
    """
}
