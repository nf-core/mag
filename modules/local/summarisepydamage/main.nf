process SUMMARISEPYDAMAGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pydamage:1.0--pyhdfd78af_0'
        : 'biocontainers/pydamage:1.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(csv)

    output:
    tuple val(meta), path("*.tsv"), emit: summary_tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    summarise_pydamage.py \\
        ${args} \\
        -i ${csv} \\
        -o ${prefix}_pydamagebins_summarised.tsv \\
        -n ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        summarisepydamage: \$(summarise_pydamage.py --version | cut -d ' ' -f 2)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    
    touch ${prefix}_pydamage_summarised.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        summarisepydamage: \$(summarise_pydamage.py --version | cut -d ' ' -f 2 )
    END_VERSIONS
    """
}
