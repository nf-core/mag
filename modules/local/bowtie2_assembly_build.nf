// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process BOWTIE2_ASSEMBLY_BUILD {
    tag "${meta.assembler}-${meta.id}"

    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.2' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/bowtie2:2.4.2--py38h1c8e9b9_1'
    } else {
        container 'quay.io/biocontainers/bowtie2:2.4.2--py38h1c8e9b9_1'
    }

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path(assembly), path('bt2_index_base*'), emit: assembly_index
    path '*.version.txt'                                    , emit: version

    script:
    def software  = getSoftwareName(task.process)
    """
    mkdir bowtie
    bowtie2-build --threads $task.cpus $assembly "bt2_index_base"
    bowtie2 --version > ${software}.version.txt
    """
}
