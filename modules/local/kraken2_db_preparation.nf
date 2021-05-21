// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process KRAKEN2_DB_PREPARATION {
    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    path db

    output:
    tuple val("${db.simpleName}"), path("database/*.k2d"), emit: db

    script:
    """
    mkdir db_tmp
    tar -xf "${db}" -C db_tmp
    mkdir database
    mv `find db_tmp/ -name "*.k2d"` database/
    """
}