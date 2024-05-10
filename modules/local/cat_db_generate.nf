process CAT_DB_GENERATE {

    conda "bioconda::cat=4.6 bioconda::diamond=2.0.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-75e2a26f10cbf3629edf2d1600db3fed5ebe6e04:eae321284604f7dabbdf121e3070bda907b91266-0' :
        'biocontainers/mulled-v2-75e2a26f10cbf3629edf2d1600db3fed5ebe6e04:eae321284604f7dabbdf121e3070bda907b91266-0' }"

    output:
    tuple env(DB_NAME), path("database/*"), path("taxonomy/*"), emit: db
    path("CAT_prepare_*.tar.gz"), optional:true               , emit: db_tar_gz
    path "versions.yml"                                       , emit: versions

    script:
    def save_db = params.save_cat_db ? "Y" : "N"
    """
    CAT prepare --fresh

    # get name/date of generated datase
    out=(*_taxonomy/)
    [[ \$out =~ (.*)_taxonomy/ ]];
    DB_NAME="CAT_prepare_\${BASH_REMATCH[1]}"

    mv *_taxonomy taxonomy
    mv *_database database
    rm database/*.nr.gz
    if [ ${save_db} = "Y" ] ; then
        tar -cf - taxonomy database | gzip > "\${DB_NAME}".tar.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CAT: \$(CAT --version | sed "s/CAT v//; s/(.*//")
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """
}
