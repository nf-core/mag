// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process MERGE_QUAST_AND_BUSCO {

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
    path(quast_bin_sum)
    path(busco_sum)

    output:
    path("quast_and_busco_summary.tsv")
    path("quast_summary.tsv")

    script:
    """
    QUAST_BIN=\$(echo \"$quast_bin_sum\" | sed 's/[][]//g')
    IFS=', ' read -r -a quast_bin <<< \"\$QUAST_BIN\"
    
    for quast_file in \"\${quast_bin[@]}\"; do
        if ! [ -f "quast_summary.tsv" ]; then 
            cp "\${quast_file}" "quast_summary.tsv"
        else
            tail -n +2 "\${quast_file}" >> "quast_summary.tsv"
        fi
    done   
    combine_tables.py $busco_sum quast_summary.tsv >quast_and_busco_summary.tsv
    """
}
