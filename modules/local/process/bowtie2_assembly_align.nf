// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BOWTIE2_ASSEMBLY_ALIGN {
    tag "${assembly_meta.assembler}-${assembly_meta.id}-${reads_meta.id}"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:"${assembly_meta.assembler}/${assembly_meta.id}_QC") }

    conda (params.enable_conda ? "bioconda::bowtie2=2.4.2 bioconda::samtools=1.11 conda-forge::pigz=2.3.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0"
    }

    input:
    tuple val(assembly_meta), path(assembly), path(index), val(reads_meta), path(reads)

    output:
    tuple val(assembly_meta), path(assembly), path("${assembly_meta.assembler}-${assembly_meta.id}-${reads_meta.id}.bam"), path("${assembly_meta.assembler}-${assembly_meta.id}-${reads_meta.id}.bam.bai"), emit: mappings
    tuple val(assembly_meta), val(reads_meta), path("*.bowtie2.log")                                                                                                                                      , emit: log
    path '*.version.txt'                                                                                                                                                                                  , emit: version

    script:
    def software = getSoftwareName(task.process)
    def name = "${assembly_meta.assembler}-${assembly_meta.id}-${reads_meta.id}"
    def input = params.single_end ? "-U \"${reads}\"" :  "-1 \"${reads[0]}\" -2 \"${reads[1]}\""
    """
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
    bowtie2 -p "${task.cpus}" -x \$INDEX $input 2> "${name}.bowtie2.log" | \
        samtools view -@ "${task.cpus}" -bS | \
        samtools sort -@ "${task.cpus}" -o "${name}.bam"
    samtools index "${name}.bam"

    if [ ${name} = "${assembly_meta.assembler}-${assembly_meta.id}-${assembly_meta.id}" ] ; then
        mv "${name}.bowtie2.log" "${assembly_meta.assembler}-${assembly_meta.id}.bowtie2.log"
    fi

    echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//' > ${software}_assembly.version.txt
    """
}
