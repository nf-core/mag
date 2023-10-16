process SPADESHYBRID {
    tag "$meta.id"

    conda "bioconda::spades=3.15.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spades:3.15.3--h95f258a_0' :
        'biocontainers/spades:3.15.3--h95f258a_0' }"

    input:
    tuple val(meta), path(long_reads), path(short_reads)

    output:
    tuple val(meta), path("SPAdesHybrid-${meta.id}_scaffolds.fasta"), emit: assembly
    path "SPAdesHybrid-${meta.id}.log"                              , emit: log
    path "SPAdesHybrid-${meta.id}_contigs.fasta.gz"                 , emit: contigs_gz
    path "SPAdesHybrid-${meta.id}_scaffolds.fasta.gz"               , emit: assembly_gz
    path "SPAdesHybrid-${meta.id}_graph.gfa.gz"                     , emit: graph
    path "versions.yml"                                , emit: versions

    script:
    def args = task.ext.args ?: ''
    maxmem = task.memory.toGiga()
    if ( params.spadeshybrid_fix_cpus == -1 || task.cpus == params.spadeshybrid_fix_cpus )
        """
        metaspades.py \
            $args \
            --threads "${task.cpus}" \
            --memory $maxmem \
            --pe1-1 ${short_reads[0]} \
            --pe1-2 ${short_reads[1]} \
            --nanopore ${long_reads} \
            -o spades
        mv spades/assembly_graph_with_scaffolds.gfa SPAdesHybrid-${meta.id}_graph.gfa
        mv spades/scaffolds.fasta SPAdesHybrid-${meta.id}_scaffolds.fasta
        mv spades/contigs.fasta SPAdesHybrid-${meta.id}_contigs.fasta
        mv spades/spades.log SPAdesHybrid-${meta.id}.log
        gzip "SPAdesHybrid-${meta.id}_contigs.fasta"
        gzip "SPAdesHybrid-${meta.id}_graph.gfa"
        gzip -c "SPAdesHybrid-${meta.id}_scaffolds.fasta" > "SPAdesHybrid-${meta.id}_scaffolds.fasta.gz"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed 's/Python //g')
            metaspades: \$(metaspades.py --version | sed "s/SPAdes genome assembler v//; s/ \\[.*//")
        END_VERSIONS
        """
    else
        error "ERROR: '--spadeshybrid_fix_cpus' was specified, but not succesfully applied. Likely this is caused by changed process properties in a custom config file."
}
