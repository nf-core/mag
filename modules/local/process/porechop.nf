// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process PORECHOP {
    tag "$meta.id"

    conda (params.enable_conda ? "bioconda::porechop=0.2.3_seqan2.1.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/porechop:0.2.3_seqan2.1.1--py36h2d50403_3"
    } else {
        container "quay.io/biocontainers/porechop:0.2.3_seqan2.1.1--py36h2d50403_3"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_porechop.fastq")

    script:
    """
    porechop -i ${reads} -t ${task.cpus} -o ${meta.id}_porechop.fastq
    """
}
