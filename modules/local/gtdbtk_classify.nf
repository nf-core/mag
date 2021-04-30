// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process GTDBTK_CLASSIFY {
    tag "${meta.assembler}-${meta.id}"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:"${meta.assembler}/${meta.id}") }

    conda (params.enable_conda ? "conda-forge::gtdbtk=1.5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/gtdbtk:1.5.0--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/gtdbtk:1.5.0--pyhdfd78af_0"
    }

    input:
    tuple val(meta), path("bins/*")
    tuple val(db_name), path("database/*")

    output:
    path 'classify/gtdbtk.*.summary.tsv'        , emit: summary
    path 'classify/gtdbtk.*.classify.tree'      , emit: tree
    path 'classify/gtdbtk.*.markers_summary.tsv', emit: markers
    path 'classify/gtdbtk.*.msa.fasta'          , emit: msa
    path 'classify/gtdbtk.*.user_msa.fasta'     , emit: user_msa
    path 'classify/gtdbtk.*.filtered.tsv'       , emit: filtered
    path 'classify/gtdbtk.log'                  , emit: log
    path 'classify/gtdbtk.warnings.log'         , emit: warnings
    path '*.version.txt'                        , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    export GTDBTK_DATA_PATH="\${PWD}/database"
    gtdbtk classify_wf --genome_dir bins --out_dir classify -x fa --cpus ${task.cpus} --min_perc_aa 5 --min_af 0.4

    gtdbtk --version | sed "s/gtdbtk: version //; s/ Copyright.*//" > ${software}.version.txt
    """
}
