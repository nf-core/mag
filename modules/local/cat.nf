process CAT {
    tag "${meta.assembler}-${meta.binner}-${meta.id}-${db_name}"

    conda "bioconda::cat=4.6 bioconda::diamond=2.0.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-75e2a26f10cbf3629edf2d1600db3fed5ebe6e04:eae321284604f7dabbdf121e3070bda907b91266-0' :
        'biocontainers/mulled-v2-75e2a26f10cbf3629edf2d1600db3fed5ebe6e04:eae321284604f7dabbdf121e3070bda907b91266-0' }"

    input:
    tuple val(meta), path("bins/*")
    tuple val(db_name), path("database/*"), path("taxonomy/*")

    output:
    path("*.names.txt")                    , emit: tax_classification
    path("raw/*.ORF2LCA.txt.gz")           , emit: orf2lca
    path("raw/*.predicted_proteins.faa.gz"), emit: faa
    path("raw/*.predicted_proteins.gff.gz"), emit: gff
    path("raw/*.log")                      , emit: log
    path("raw/*.bin2classification.txt.gz"), emit: tax_classification_taxids
    path "versions.yml"                    , emit: versions

    script:
    def official_taxonomy = params.cat_official_taxonomy ? "--only_official" : ""
    """
    CAT bins -b "bins/" -d database/ -t taxonomy/ -n "${task.cpus}" -s .fa --top 6 -o "${meta.assembler}-${meta.binner}-${meta.id}" --I_know_what_Im_doing
    CAT add_names -i "${meta.assembler}-${meta.binner}-${meta.id}.ORF2LCA.txt" -o "${meta.assembler}-${meta.binner}-${meta.id}.ORF2LCA.names.txt" -t taxonomy/ ${official_taxonomy}
    CAT add_names -i "${meta.assembler}-${meta.binner}-${meta.id}.bin2classification.txt" -o "${meta.assembler}-${meta.binner}-${meta.id}.bin2classification.names.txt" -t taxonomy/ ${official_taxonomy}

    mkdir raw
    mv *.ORF2LCA.txt *.predicted_proteins.faa *.predicted_proteins.gff *.bin2classification.txt *.log raw/
    cp *.bin2classification.names.txt raw/

    gzip "raw/${meta.assembler}-${meta.binner}-${meta.id}.ORF2LCA.txt" \
        "raw/${meta.assembler}-${meta.binner}-${meta.id}.concatenated.predicted_proteins.faa" \
        "raw/${meta.assembler}-${meta.binner}-${meta.id}.concatenated.predicted_proteins.gff" \
        "raw/${meta.assembler}-${meta.binner}-${meta.id}.bin2classification.txt" \
        "${meta.assembler}-${meta.binner}-${meta.id}.ORF2LCA.names.txt" \

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CAT: \$(CAT --version | sed "s/CAT v//; s/(.*//")
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """
}
