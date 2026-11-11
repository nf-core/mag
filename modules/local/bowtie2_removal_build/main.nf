process BOWTIE2_REMOVAL_BUILD {
    tag "${fasta}"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bowtie2:2.4.2--py38h1c8e9b9_1'
        : 'biocontainers/bowtie2:2.4.2--py38h1c8e9b9_1'}"

    input:
    path fasta

    output:
    path "*.{bt2,bt2l}", emit: index
    path "versions.yml", emit: versions

    script:
    """
    mkdir bowtie
    bowtie2-build --threads ${task.cpus} ${fasta} ${fasta.simpleName}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """
}
