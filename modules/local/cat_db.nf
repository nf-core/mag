// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process CAT_DB {
    tag "${database.baseName}"

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    path(database)

    output:
    tuple val("${database.toString().replace(".tar.gz", "")}"), path("database/*"), path("taxonomy/*")

    script:
    """
    mkdir catDB
    tar -xf ${database} -C catDB
    mv `find catDB/ -type d -name "*taxonomy*"` taxonomy/
    mv `find catDB/ -type d -name "*CAT_database*"` database/
    """
}
