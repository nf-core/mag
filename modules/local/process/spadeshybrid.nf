// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SPADESHYBRID {
    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::spades=3.13.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/spades:3.13.1--h2d02072_2"
    } else {
        container "quay.io/biocontainers/spades:3.13.1--h2d02072_2"
    }

    input:
    tuple val(meta), path(long_reads), path(short_reads)

    output:
    tuple val("SPAdesHybrid"), val(meta), path("${meta.id}_scaffolds.fasta"), emit: assembly
    file("${meta.id}.log")
    file("${meta.id}_contigs.fasta.gz")
    file("${meta.id}_scaffolds.fasta.gz")
    file("${meta.id}_graph.gfa.gz")

    script:
    maxmem = task.memory.toGiga()
    if ( !params.spadeshybrid_fix_cpus || task.cpus == params.spadeshybrid_fix_cpus )
        """
        metaspades.py \
            --threads "${task.cpus}" \
            --memory $maxmem \
            --pe1-1 ${short_reads[0]} \
            --pe1-2 ${short_reads[1]} \
            --nanopore ${long_reads} \
            -o spades
        mv spades/assembly_graph_with_scaffolds.gfa ${meta.id}_graph.gfa
        mv spades/scaffolds.fasta ${meta.id}_scaffolds.fasta
        mv spades/contigs.fasta ${meta.id}_contigs.fasta
        mv spades/spades.log ${meta.id}.log
        gzip "${meta.id}_contigs.fasta"
        gzip "${meta.id}_graph.gfa"
        gzip -c "${meta.id}_scaffolds.fasta" > "${meta.id}_scaffolds.fasta.gz"
        """
    else
        error "ERROR: '--spadeshybrid_fix_cpus' was specified, but not succesfully applied. Likely this is caused by changed process properties in a custom config file."
}
