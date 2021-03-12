// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process KRAKEN2_DB_PREPARATION {
    input:
    path db

    output:
    tuple val("${db.baseName}"), path("*/*.k2d")

    script:
    """
    tar -xf "${db}"
    """
}