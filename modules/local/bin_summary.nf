// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process BIN_SUMMARY {

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
    path(bin_depths)
    path(busco_sum)
    path(quast_sum)
    path(gtdbtk_sum)

    output:
    path("bin_summary.tsv"), emit: summary

    script:
    def busco_summary  = busco_sum.sort().size() > 0 ?  "--busco_summary ${busco_sum}" : ""
    def quast_summary  = quast_sum.sort().size() > 0 ?  "--quast_summary ${quast_sum}" : ""
    def gtdbtk_summary = gtdbtk_sum.sort().size() > 0 ? "--gtdbtk_summary ${gtdbtk_sum}" : ""
    """
    combine_tables.py --depths_summary ${bin_depths} \
                      $busco_summary \
                      $quast_summary \
                      $gtdbtk_summary \
                      --out bin_summary.tsv
    """
}
