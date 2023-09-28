process KRAKEN2_DB_PREPARATION {

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'docker.io/ubuntu:20.04' }"

    input:
    path db

    output:
    tuple val("${db.simpleName}"), path("database/*.k2d"), emit: db
    path "versions.yml"                                  , emit: versions

    script:
    """
    if [[ -d ${db} ]]; then
        if [[ ${db} != database ]]; then
            ln -sr ${db} database
        fi

        # Make sure {hash,opts,taxo}.k2d are found in directory input
        if [[ \$(find database/ -name "*.k2d" | wc -l) -lt 3 ]]; then
            error "ERROR: Kraken2 requires '{hash,opts,taxo}.k2d' files."
        fi
    else
        mkdir db_tmp
        tar -xf "${db}" -C db_tmp
        mkdir database
        mv `find db_tmp/ -name "*.k2d"` database/
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(tar --version 2>&1 | sed -n 1p | sed 's/tar (GNU tar) //')
    END_VERSIONS
    """
}
