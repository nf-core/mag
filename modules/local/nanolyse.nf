// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process NANOLYSE {
    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::nanolyse=1.1.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/nanolyse:1.1.0--py36_1"
    } else {
        container "quay.io/biocontainers/nanolyse:1.1.0--py36_1"
    }

    input:
    tuple val(meta), path(reads)
    path nanolyse_db

    output:
    tuple val(meta), path("${meta.id}_nanolyse.fastq.gz"), emit: reads
    path  "${meta.id}_nanolyse.log"                      , emit: log
    path '*.version.txt'                                 , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    cat ${reads} | NanoLyse --reference $nanolyse_db | gzip > ${meta.id}_nanolyse.fastq.gz
    echo "NanoLyse reference: $params.lambda_reference" >${meta.id}_nanolyse.log
    cat ${reads} | echo "total reads before NanoLyse: \$((`wc -l`/4))" >>${meta.id}_nanolyse.log
    gunzip -c ${meta.id}_nanolyse.fastq.gz | echo "total reads after NanoLyse: \$((`wc -l`/4))" >> ${meta.id}_nanolyse.log

    NanoLyse --version | sed -e "s/NanoLyse //g" > ${software}.version.txt
    """
}
