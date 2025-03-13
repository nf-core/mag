process BUSCO {
    tag "${bin}"

    conda "bioconda::busco=5.4.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.4.3--pyhdfd78af_0':
        'biocontainers/busco:5.4.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bin)
    tuple val(db_meta), path(db)

    output:
    tuple val(meta), path("short_summary.domain.*.${bin}.txt")          , optional:true , emit: summary_domain
    tuple val(meta), path("short_summary.specific_lineage.*.${bin}.txt"), optional:true , emit: summary_specific
    tuple env('MOST_SPEC_DB'), path('busco_downloads/')                  , optional:true , emit: busco_downloads
    path("${bin}_busco.log")
    path("${bin}_busco.err")
    path("${bin}_buscos.*.faa.gz")                                      , optional:true
    path("${bin}_buscos.*.fna.gz")                                      , optional:true
    path("${bin}_prodigal.gff")                                         , optional:true , emit: prodigal_genes
    tuple val(meta), path("${bin}_busco.failed_bin.txt")                , optional:true , emit: failed_bin
    path "versions.yml"                                                                 , emit: versions

    script:
    def args = task.ext.args ?: ''
    def cp_augustus_config = workflow.profile.toString().indexOf("conda") != -1 ? "N" : "Y"
    def lineage_dataset_provided = "${db_meta.lineage}"
    def busco_clean = params.busco_clean ? "Y" : "N"

    def p = params.busco_auto_lineage_prok ? "--auto-lineage-prok" : "--auto-lineage"
    if ( "${lineage_dataset_provided}" == "Y" ) {
        p = "--lineage_dataset dataset/${db}"
    } else if ( "${lineage_dataset_provided}" == "N" ) {
        p += " --offline --download_path ${db}"
    } else {
        lineage_dataset_provided = ""
    }
    """
    run_busco.sh \\
        "${p}" \\
        "${cp_augustus_config}" \\
        "${db}" \\
        "${bin}" \\
        ${task.cpus} \\
        "${lineage_dataset_provided}" \\
        "${busco_clean}" \\
        "${args}"

    MOST_SPEC_DB=\$(<info_most_spec_db.txt)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        R: \$(R --version 2>&1 | sed -n 1p | sed 's/R version //' | sed 's/ (.*//')
        busco: \$(busco --version 2>&1 | sed 's/BUSCO //g')
    END_VERSIONS
    """
}
