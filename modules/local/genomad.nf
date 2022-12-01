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
    // TODO: Add additional output channels here
    tuple val(meta), path("${meta.id}_summary/*.tsv"), emit: tsv
    path("*.log")                                    , emit: log
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def clean_db= db.toString() - ".tar.gz"
    def decompress_db= db.toString() == clean_db ? "" : "mkdir $clean_db && tar xf $db -C $clean_db --strip-components 1"
    """
    $decompress_db

    genomad end-to-end $args $fasta . $clean_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(echo \$(genomad --version 2>&1) | sed 's/^.*genomad, version //')
    END_VERSIONS
    """
}
