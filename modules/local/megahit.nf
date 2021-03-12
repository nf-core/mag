// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MEGAHIT {
    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::megahit=1.2.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/megahit:1.2.7--h8b12597_0"
    } else {
        container "quay.io/biocontainers/megahit:1.2.7--h8b12597_0"
    }

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("MEGAHIT/${meta.id}.contigs.fa"), emit: assembly
    path "MEGAHIT/*.log"
    path "MEGAHIT/${meta.id}.contigs.fa.gz"
    path '*.version.txt'                                  , emit: version

    script:
    def software = getSoftwareName(task.process)
    def input = params.single_end ? "-r \"" + reads1.join(",") + "\"" : "-1 \"" + reads1.join(",") + "\" -2 \"" + reads2.join(",") + "\""
    mem = task.memory.toBytes()
    if ( !params.megahit_fix_cpu_1 || task.cpus == 1 )
        """
        megahit ${params.megahit_options} -t "${task.cpus}" -m $mem $input -o MEGAHIT --out-prefix "${meta.id}"
        gzip -c "MEGAHIT/${meta.id}.contigs.fa" > "MEGAHIT/${meta.id}.contigs.fa.gz"

        megahit --version | sed "s/MEGAHIT v//" > ${software}.version.txt
        """
    else
        error "ERROR: '--megahit_fix_cpu_1' was specified, but not succesfully applied. Likely this is caused by changed process properties in a custom config file."
}
