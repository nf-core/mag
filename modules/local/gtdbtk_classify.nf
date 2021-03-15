// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process GTDBTK_CLASSIFY {
    tag "${meta.assembler}-${meta.id}"

    conda (params.enable_conda ? "conda-forge::gtdbtk=1.4.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gtdbtk:1.4.1--py_1"
    } else {
        container "quay.io/biocontainers/gtdbtk:1.4.1--py_1"
    }

    input:
    tuple val(meta), path("bins/*")
    tuple val(db_name), path("database/*")

    output:
    tuple val(meta), path("classify/*"), emit: taxonomy
    path '*.version.txt'               , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    export GTDBTK_DATA_PATH="\${PWD}/database"
    gtdbtk classify_wf --genome_dir bins --out_dir classify -x fa --cpus ${task.cpus} --min_perc_aa 5 --min_af 0.4

    gtdbtk --version | sed "s/gtdbtk: version //; s/ Copyright.*//" > ${software}.version.txt
    """
}
