// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process BOWTIE2_REMOVAL_BUILD {
    tag "$fasta"

    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.2' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/bowtie2:2.4.2--py38h1c8e9b9_1'
    } else {
        container 'quay.io/biocontainers/bowtie2:2.4.2--py38h1c8e9b9_1'
    }

    input:
    path fasta

    output:
    path 'bt2_index_base*', emit: index
    path '*.version.txt'  , emit: version

    script:
    def software  = getSoftwareName(task.process)
    """
    mkdir bowtie
    bowtie2-build --threads $task.cpus $fasta "bt2_index_base"
    bowtie2 --version > ${software}_removal.version.txt
    """
}
