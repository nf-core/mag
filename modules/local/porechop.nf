process PORECHOP {
    tag "$meta.id"

    conda "bioconda::porechop=0.2.3_seqan2.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/porechop:0.2.3_seqan2.1.1--py36h2d50403_3' :
        'quay.io/biocontainers/porechop:0.2.3_seqan2.1.1--py36h2d50403_3' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_porechop.fastq")  , emit: reads
    path "versions.yml"                                 , emit: versions

    script:
    """
    porechop -i ${reads} -t ${task.cpus} -o ${meta.id}_porechop.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        porechop: \$(porechop --version)
    END_VERSIONS
    """
}
