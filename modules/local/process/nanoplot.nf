// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process NANOPLOT {
    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::nanoplot=1.26.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/nanoplot:1.26.3--py_0"
    } else {
        container "quay.io/biocontainers/nanoplot:1.26.3--py_0"
    }

    input:
    tuple val(meta), path(reads)

    output:
    file '*.png'
    file '*.html'
    file '*.txt'

    script:
    def prefix = options.suffix ? "-p ${options.suffix}_" : ''
    def title  = options.suffix ? "${meta.id}_${options.suffix}" : "${meta.id}"
    """
    NanoPlot -t ${task.cpus} \
             ${prefix} \
             --title ${title} \
             -c darkblue \
             --fastq ${reads}
    """
}