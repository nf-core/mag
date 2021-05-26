// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process QUAST_BINS_SUMMARY {

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    path(summaries)

    output:
    path("quast_summary.tsv"), emit: summary

    script:
    """
    QUAST_BIN=\$(echo \"$summaries\" | sed 's/[][]//g')
    IFS=', ' read -r -a quast_bin <<< \"\$QUAST_BIN\"
    for quast_file in \"\${quast_bin[@]}\"; do
        if ! [ -f "quast_summary.tsv" ]; then 
            cp "\${quast_file}" "quast_summary.tsv"
        else
            tail -n +2 "\${quast_file}" >> "quast_summary.tsv"
        fi
    done
    """
}
