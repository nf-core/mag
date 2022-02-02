process POOL_PAIRED_READS {
    tag "$meta.id"

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("pooled_${meta.id}_*.fastq.gz"), emit: reads
    path "versions.yml"                                  , emit: versions

    script:
    """
    cat ${reads1} > "pooled_${meta.id}_1.fastq.gz"
    cat ${reads2} > "pooled_${meta.id}_2.fastq.gz"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """
}
