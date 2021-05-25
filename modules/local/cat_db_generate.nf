// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process CAT_DB_GENERATE {

    publishDir "${params.outdir}",
        mode: 'move',
        saveAs: { filename -> params.save_cat_db ? saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') : null }

    conda (params.enable_conda ? "bioconda::cat=4.6 bioconda::diamond=2.0.6" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-75e2a26f10cbf3629edf2d1600db3fed5ebe6e04:eae321284604f7dabbdf121e3070bda907b91266-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-75e2a26f10cbf3629edf2d1600db3fed5ebe6e04:eae321284604f7dabbdf121e3070bda907b91266-0"
    }

    output:
    tuple env(DB_NAME), path("database/*"), path("taxonomy/*"), emit: db
    path("CAT_prepare_*.tar.gz"), optional:true               , emit: db_tar_gz
    path '*.version.txt'                                      , emit: version

    script:
    def software = getSoftwareName(task.process)
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

    CAT --version | sed "s/CAT v//; s/(.*//" > ${software}.version.txt
    """
}
