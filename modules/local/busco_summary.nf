// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process BUSCO_SUMMARY {

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "conda-forge::python:3.6.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.6.7"
    } else {
        container "quay.io/biocontainers/python:3.6.7"
    }

    input:
    path "short_summary.*.txt"

    output:
    path "busco_summary.txt"

    script:
    """
    summary_busco.py short_summary.*.txt > busco_summary.txt
    """
}
