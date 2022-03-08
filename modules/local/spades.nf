process SPADES {
    tag "$meta.id"

    conda (params.enable_conda ? "bioconda::spades=3.15.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spades:3.15.3--h95f258a_0' :
        'quay.io/biocontainers/spades:3.15.3--h95f258a_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("SPAdes-${meta.id}_scaffolds.fasta"), emit: assembly
    path "SPAdes-${meta.id}.log"                              , emit: log
    path "SPAdes-${meta.id}_contigs.fasta.gz"                 , emit: contigs_gz
    path "SPAdes-${meta.id}_scaffolds.fasta.gz"               , emit: assembly_gz
    path "SPAdes-${meta.id}_graph.gfa.gz"                     , emit: graph
    path "versions.yml"                                , emit: versions

    script:
    def args = task.ext.args ?: ''
    maxmem = task.memory.toGiga()
    if ( params.spades_fix_cpus == -1 || task.cpus == params.spades_fix_cpus )
        """
        metaspades.py \
            $args \
            --threads "${task.cpus}" \
            --memory $maxmem \
            --pe1-1 ${reads[0]} \
            --pe1-2 ${reads[1]} \
            -o spades
        mv spades/assembly_graph_with_scaffolds.gfa SPAdes-${meta.id}_graph.gfa
        mv spades/scaffolds.fasta SPAdes-${meta.id}_scaffolds.fasta
        mv spades/contigs.fasta SPAdes-${meta.id}_contigs.fasta
        mv spades/spades.log SPAdes-${meta.id}.log
        gzip "SPAdes-${meta.id}_contigs.fasta"
        gzip "SPAdes-${meta.id}_graph.gfa"
        gzip -c "SPAdes-${meta.id}_scaffolds.fasta" > "SPAdes-${meta.id}_scaffolds.fasta.gz"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version 2>&1 | sed 's/Python //g')
            metaspades: \$(metaspades.py --version | sed "s/SPAdes genome assembler v//; s/ \\[.*//")
        END_VERSIONS
        """
    else
        error "ERROR: '--spades_fix_cpus' was specified, but not succesfully applied. Likely this is caused by changed process properties in a custom config file."
}
