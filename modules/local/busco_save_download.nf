// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process BUSCO_SAVE_DOWNLOAD {
    // execute sequentially to avoid artefacts when saving files for multiple busco instances
    maxForks 1

    // do not overwrite existing files which were already saved for other busco runs
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        overwrite: false,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    path(busco_downloads)

    output:
    path('busco_downloads/**', includeInputs: true)

    script:
    """
    """
}
