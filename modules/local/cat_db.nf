process CAT_DB {
    tag "${database.baseName}"

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    path(database)

    output:
    tuple val("${database.toString().replace(".tar.gz", "")}"), path("database/*"), path("taxonomy/*"), emit: db
    path "versions.yml"                                                                               , emit: versions

    script:
    """
    mkdir catDB
    tar -xf ${database} -C catDB
    mv `find catDB/ -type d -name "*taxonomy*"` taxonomy/
    mv `find catDB/ -type d -name "*database*"` database/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(tar --version 2>&1 | sed -n 1p | sed 's/tar (GNU tar) //')
    END_VERSIONS
    """
}
