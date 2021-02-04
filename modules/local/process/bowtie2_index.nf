// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BOWTIE2_INDEX {
    tag "$fasta"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? 'bioconda::bowtie2=2.4.2' : null) // TODO use previous version, update tools separately!!!
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
    bowtie2 --version > ${software}.version.txt
    """
}
