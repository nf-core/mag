process GENOMAD {
    tag "${meta.id}-${db.simpleName}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::genomad=1.2.0" : null)
    // TODO: Replace the container below with the followng once the container is available
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/genomad=1.2.0':
    //     'quay.io/biocontainers/genomad=1.2.0' }"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://philpalmer/genomad:1.2.0':
        'philpalmer/genomad:1.2.0' }"

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    path("${prefix}_annotate/${prefix}_taxonomy.tsv")                                  , emit: taxonomy_tsv
    path("${prefix}_aggregated_classification/${prefix}_aggregated_classification.tsv"), emit: aggregated_classification_tsv
    path("${prefix}_summary/${prefix}_virus_summary.tsv")                              , emit: virus_summary_tsv
    path("${prefix}_summary/${prefix}_plasmid_summary.tsv")                            , emit: plasmid_summary_tsv
    path("${prefix}_summary/${prefix}_virus_genes.tsv")                                , emit: viruses_genes_tsv
    path("${prefix}_summary/${prefix}_plasmid_genes.tsv")                              , emit: plasmids_genes_tsv
    path("${prefix}_summary/${prefix}_virus.fna")                                      , emit: viruses_fna
    path("${prefix}_summary/${prefix}_plasmid.fna")                                    , emit: plasmids_fna
    path("${prefix}_summary/${prefix}_virus_proteins.faa")                             , emit: viruses_faa
    path("${prefix}_summary/${prefix}_plasmid_proteins.faa")                           , emit: plasmids_faa
    path("${prefix}.log")                                                              , emit: log
    path "versions.yml"                                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def clean_db = db.toString() - ".tar.gz"
    def decompress_db = db.toString() == clean_db ? "" : "rm -rf $clean_db && mkdir $clean_db && tar xf $db -C $clean_db --strip-components 1"
    prefix = fasta.toString().replaceAll(/\.gz$/, '').replaceAll(/(\.fasta|\.fa)$/, '')
    """
    $decompress_db

    genomad end-to-end $args $fasta . $clean_db

    ls -1rt ${prefix}*.log | head -n6 | xargs cat > ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(echo \$(genomad --version 2>&1) | sed 's/^.*genomad, version //')
    END_VERSIONS
    """
}
