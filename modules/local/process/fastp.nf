// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

/*
 * FASTP
 */
process FASTP {
    tag "$meta.id"
    publishDir "${params.outdir}/",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::fastp=0.20.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0"
    } else {
        container "quay.io/biocontainers/fastp:0.20.1--h8b12597_0"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_trimmed*.fastq.gz"), emit: reads
    path  "fastp.html"
    path  "fastp.json"
    path  "*.version.txt"                                , emit: version

    script:
    def software = getSoftwareName(task.process)
    def pe_input = meta.single_end ? '' :  "-I \"${reads[1]}\""
    def pe_output1 = meta.single_end ? "-o \"${meta.id}_trimmed.fastq.gz\"" :  "-o \"${meta.id}_trimmed_R1.fastq.gz\""
    def pe_output2 = meta.single_end ? '' :  "-O \"${meta.id}_trimmed_R2.fastq.gz\""
    """
    fastp -w ${task.cpus} \
          -q ${params.mean_quality} \
          --cut_by_quality5 \
          --cut_by_quality3 \
          --cut_mean_quality ${params.trimming_quality} \
          --adapter_sequence=${params.adapter_forward} \
          --adapter_sequence_r2=${params.adapter_reverse} \
          -i "${reads[0]}" $pe_input \
          $pe_output1 $pe_output2

    echo \$(fastp --version 2>&1) | sed 's/fastp //' > ${software}.version.txt
    """
}
