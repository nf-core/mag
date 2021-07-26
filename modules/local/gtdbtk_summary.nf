// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process GTDBTK_SUMMARY {

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::pandas=1.1.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pandas:1.1.5"
    } else {
        container "quay.io/biocontainers/pandas:1.1.5"
    }

    input:
    path(qc_discarded_bins)
    path(gtdbtk_summaries)
    path(filtered_bins)
    path(failed_bins)

    output:
    path "gtdbtk_summary.tsv", emit: summary

    script:
    def discarded = qc_discarded_bins.sort().size() > 0 ? "--qc_discarded_bins ${qc_discarded_bins}" : ""
    def summaries = gtdbtk_summaries.sort().size() > 0 ?  "--summaries ${gtdbtk_summaries}" : ""
    def filtered  = filtered_bins.sort().size() > 0 ?     "--filtered_bins ${filtered_bins}" : ""
    def failed    = failed_bins.sort().size() > 0 ?       "--failed_bins ${failed_bins}" : ""
    """
    summary_gtdbtk.py $options.args $discarded $summaries $filtered $failed --out gtdbtk_summary.tsv
    """
}
