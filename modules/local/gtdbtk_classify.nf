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
    path "classify/gtdbtk.${meta.assembler}-${meta.id}.*.summary.tsv"        , emit: summary
    path "classify/gtdbtk.${meta.assembler}-${meta.id}.*.classify.tree"      , emit: tree
    path "classify/gtdbtk.${meta.assembler}-${meta.id}.*.markers_summary.tsv", emit: markers
    path "classify/gtdbtk.${meta.assembler}-${meta.id}.*.msa.fasta"          , emit: msa
    path "classify/gtdbtk.${meta.assembler}-${meta.id}.*.user_msa.fasta"     , emit: user_msa
    path "classify/gtdbtk.${meta.assembler}-${meta.id}.*.filtered.tsv"       , emit: filtered
    path "classify/gtdbtk.${meta.assembler}-${meta.id}.log"                  , emit: log
    path "classify/gtdbtk.${meta.assembler}-${meta.id}.warnings.log"         , emit: warnings
    path "classify/gtdbtk.${meta.assembler}-${meta.id}.failed_genomes.tsv"   , emit: failed
    path '*.version.txt'                                                     , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    export GTDBTK_DATA_PATH="\${PWD}/database"
    gtdbtk classify_wf --genome_dir bins --prefix "gtdbtk.${meta.assembler}-${meta.id}" --out_dir classify -x fa --cpus ${task.cpus} --min_perc_aa ${params.gtdbtk_min_perc_aa} --min_af ${params.gtdbtk_min_af}

    mv classify/gtdbtk.log "classify/gtdbtk.${meta.assembler}-${meta.id}.log"
    mv classify/gtdbtk.warnings.log "classify/gtdbtk.${meta.assembler}-${meta.id}.warnings.log"
    gtdbtk --version | sed "s/gtdbtk: version //; s/ Copyright.*//" > ${software}.version.txt
    """
}
