// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process BUSCO_SUMMARY {

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "conda-forge::pandas=1.1.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pandas:1.1.5"
    } else {
        container "quay.io/biocontainers/pandas:1.1.5"
    }

    input:
    path(summaries_domain)
    path(summaries_specific)
    path(failed_bins)

    output:
    path "busco_summary.tsv", emit: summary

    script:
    def auto = params.busco_reference ? "" : "-a"
    def ss = summaries_specific.sort().size() > 0 ? "-ss ${summaries_specific}" : ""
    def sd = summaries_domain.sort().size() > 0 ? "-sd ${summaries_domain}" : ""
    def f = ""
    if (!params.busco_reference && failed_bins.sort().size() > 0)
        f = "-f ${failed_bins}"
    """
    summary_busco.py $auto $ss $sd $f -o busco_summary.tsv
    """
}

