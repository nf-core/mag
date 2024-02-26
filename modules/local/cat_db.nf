process CAT_DB {
    tag "${database.baseName}"

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path(database)

    output:
    tuple val("${database.toString().replace(".tar.gz", "")}"), path("database/*"), path("taxonomy/*"), emit: db
    path "versions.yml"                                                                               , emit: versions

    script:
    """
    if [[ ${database} != *.tar.gz ]]; then
        ln -sr `find ${database}/ -type d -name "*taxonomy*"` taxonomy
        ln -sr `find ${database}/ -type d -name "*database*"` database
    else
        mkdir catDB
        tar -xf ${database} -C catDB
        mv `find catDB/ -type d -name "*taxonomy*"` taxonomy/
        mv `find catDB/ -type d -name "*database*"` database/
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(tar --version 2>&1 | sed -n 1p | sed 's/tar (GNU tar) //')
    END_VERSIONS
    """
}
