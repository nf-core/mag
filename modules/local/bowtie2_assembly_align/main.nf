process BOWTIE2_ASSEMBLY_ALIGN {
    tag "${assembly_meta.assembler}-${assembly_meta.id}-${reads_meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0' :
        'biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:577a697be67b5ae9b16f637fd723b8263a3898b3-0' }"

    input:
    tuple val(assembly_meta), path(assembly), path(index), val(reads_meta), path(reads)

    output:
    tuple val(assembly_meta), path(assembly), path("${assembly_meta.assembler}-${assembly_meta.id}-${reads_meta.id}.bam"), path("${assembly_meta.assembler}-${assembly_meta.id}-${reads_meta.id}.bam.bai"), emit: mappings
    tuple val(assembly_meta), val(reads_meta), path("*.bowtie2.log")                                                                                                                                      , emit: log
    path "versions.yml"                                                                                                                                                                                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    def name = "${assembly_meta.assembler}-${assembly_meta.id}-${reads_meta.id}"
    def input = params.single_end ? "-U \"${reads}\"" :  "-1 \"${reads[0]}\" -2 \"${reads[1]}\""
    """
    INDEX=`find -L ./ -name "*.rev.1.bt2l" -o -name "*.rev.1.bt2" | sed 's/.rev.1.bt2l//' | sed 's/.rev.1.bt2//'`
    bowtie2 \\
        -p "${task.cpus}" \\
        -x \$INDEX \\
        $args \\
        $input \\
        2> "${name}.bowtie2.log" | \
        samtools view -@ "${task.cpus}" -bS | \
        samtools sort -@ "${task.cpus}" -o "${name}.bam"
    samtools index "${name}.bam"

    if [ ${name} = "${assembly_meta.assembler}-${assembly_meta.id}-${assembly_meta.id}" ] ; then
        mv "${name}.bowtie2.log" "${assembly_meta.assembler}-${assembly_meta.id}.bowtie2.log"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
