process POOL_SINGLE_READS {
    tag "$meta.id"

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("pooled_${meta.id}.fastq.gz"), emit: reads
    path "versions.yml"                                , emit: versions

    script:
    """
    cat ${reads} > "pooled_${meta.id}.fastq.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """
}
