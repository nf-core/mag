process GTDBTK_DB_PREPARATION {
    tag "${database}"

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path(database)

    output:
    tuple val("${database.toString().replace(".tar.gz", "")}"), path("database/*"), emit: db

    script:
    """
    mkdir database
    tar -xzf ${database} -C database --strip 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tar: \$(tar --version 2>&1 | sed -n 1p | sed 's/tar (GNU tar) //')
    END_VERSIONS
    """
}

/* Mount the database image and inspect directory structure.  Generate a container mount directive for each
 * depth-2 file and directory.  This is for compatability with GTDBTK_DB_PREPARATION which removes first path
 * component during tar extraction.
 */
process GTDBTK_DB_INSPECT {
    tag "${database}"

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    // Variable NXF_TASK_WORKDIR will be evaluated in container's environment.
    containerOptions "-B $database_image:\$NXF_TASK_WORKDIR/database:image-src=/"

    input:
    val(database_image)

    output:
    stdout emit: mount

    script:
    """
    find -L database -mindepth 2 -maxdepth 2 \
    -printf \"-B $database_image:\\\$NXF_TASK_WORKDIR/database/%f:image-src=/%P \"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        find: \$(find --version 2>&1 | sed -n 1p | sed 's/find (GNU findutils) //')
    END_VERSIONS
    """
}
