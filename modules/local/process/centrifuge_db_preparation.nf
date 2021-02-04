// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process CENTRIFUGE_DB_PREPARATION {
    input:
    path db

    output:
    tuple val("${db.toString().replace(".tar.gz", "")}"), path("*.cf")

    script:
    """
    tar -xf "${db}"
    """
}
