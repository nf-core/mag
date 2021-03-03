// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BOWTIE2_ASSEMBLY {
    tag "${assembler}-${assembly_id}"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:"${assembler}/${assembly_id}_QC") }

    conda (params.enable_conda ? "bioconda::bowtie2=2.4.2=py38h1c8e9b9_1 bioconda::samtools=1.11=h6270b1f_0 conda-forge::pigz=2.3.4=hed695b0_1" : null)  // TODO ?
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0"
    }

    input:
    tuple val(assembler), val(assembly_id), path(assembly), path(index), val(reads_id), path(reads)

    output:
    tuple val(assembler), val(assembly_id), path(assembly), path("${assembler}-${assembly_id}-${reads_id}.bam"), path("${assembler}-${assembly_id}-${reads_id}.bam.bai"), emit: mappings
    tuple val(assembler), val(assembly_id), val(reads_id), path("*.bowtie2.log")                                                                                        , emit: log
    path '*.version.txt'                                                                                                                                                , emit: version

    script:
    def software = getSoftwareName(task.process)
    def name = "${assembler}-${assembly_id}-${reads_id}"
    def input = params.single_end ? "-U \"${reads}\"" :  "-1 \"${reads[0]}\" -2 \"${reads[1]}\""
    """
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
    bowtie2 -p "${task.cpus}" -x \$INDEX $input 2> "${name}.bowtie2.log" | \
        samtools view -@ "${task.cpus}" -bS | \
        samtools sort -@ "${task.cpus}" -o "${name}.bam"
    samtools index "${name}.bam"

    if [ ${name} = "${assembler}-${assembly_id}-${assembly_id}" ] ; then
        mv "${name}.bowtie2.log" "${assembler}-${assembly_id}.bowtie2.log"
    fi

    echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//' > ${software}_assembly.version.txt
    """
}
