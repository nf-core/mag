// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BUSCO_DB_PREPARATION {
    tag "${database.baseName}"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> params.save_busco_reference ? saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') : null }

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    path(database)

    output:
    path "buscodb/*", emit: db
    path database

    script:
    """
    mkdir buscodb
    tar -xf ${database} -C buscodb
    """
}
