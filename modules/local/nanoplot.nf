// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

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
    path '*.png'
    path '*.html'
    path '*.txt'
    path '*.version.txt', emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix = options.suffix ? "-p ${options.suffix}_" : ''
    def title  = options.suffix ? "${meta.id}_${options.suffix}" : "${meta.id}"
    """
    NanoPlot -t ${task.cpus} \
             ${prefix} \
             --title ${title} \
             -c darkblue \
             --fastq ${reads}
    NanoPlot --version | sed -e "s/NanoPlot //g" > ${software}.version.txt
    """
}