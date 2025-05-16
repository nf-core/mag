process CAT_DB {
    tag "${database.baseName}"

    conda "conda-forge::sed=4.7"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/ubuntu:20.04'
        : 'nf-core/ubuntu:20.04'}"

    input:
    path database

    output:
    tuple val(name), path("database/"), emit: db
    tuple val(name), path("taxonomy/"), emit: taxonomy
    path "versions.yml", emit: versions

    script:
    name = "${database.toString().replace(".tar.gz", "")}"
    """
    if [[ ${database} != *.tar.gz ]]; then
        ln -sr ${database}/tax taxonomy
        ln -sr ${database}/db database
    else
        mkdir catDB
        tar -xf ${database} -C catDB
        mv `find catDB/ -type d -name "*tax*"` taxonomy/
        mv `find catDB/ -type d -name "*db*"` database/
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(tar --version 2>&1 | sed -n 1p | sed 's/tar (GNU tar) //')
    END_VERSIONS
    """
}
