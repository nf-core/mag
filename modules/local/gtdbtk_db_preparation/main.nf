process GTDBTK_DB_PREPARATION {
    tag "${database}"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c2/c262fc09eca59edb5a724080eeceb00fb06396f510aefb229c2d2c6897e63975/data'
        : 'community.wave.seqera.io/library/coreutils:9.5--ae99c88a9b28c264'}"

    input:
    path database

    output:
    tuple val("${database.toString().replace(".tar.gz", "")}"), path("database/*"), emit: db
    path "versions.yml", emit: versions

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
