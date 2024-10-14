process CONVERT_DEPTHS {
    label 'process_single'
    tag "${meta.id}"
    conda "bioconda::bioawk=1.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bioawk:1.0--hed695b0_5'
        : 'biocontainers/bioawk:1.0--hed695b0_5'}"

    input:
    tuple val(meta), path(fasta), path(depth)

    output:
    tuple val(meta), path(fasta), val([]), path("*_mb2_depth.txt"), emit: output
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gunzip -f ${depth}
    bioawk -t '{ { if (NR > 1) { { print \$1, \$3 } } } }' ${depth.toString() - '.gz'} > ${prefix}_mb2_depth.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioawk: \$(bioawk --version | cut -f 3 -d ' ' )
    END_VERSIONS
    """
}
