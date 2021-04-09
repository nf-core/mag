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
    path(summaries)
    path(failed_bins)

    output:
    path "busco_summary.txt", emit: summary

    script:
    def s = summaries.size() > 0 ? "-s ${summaries}" : ""
    def f = ""
    if (!params.busco_reference && failed_bins.size() > 0)
        f = "-f ${failed_bins}"
    """
    summary_busco.py $s $f > busco_summary.txt
    """
}

