process FILTLONG {
    tag "$meta.id"

    conda "bioconda::filtlong=0.2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/filtlong:0.2.0--he513fc3_3' :
        'quay.io/biocontainers/filtlong:0.2.0--he513fc3_3' }"

    input:
    tuple val(meta), path(long_reads), path(short_reads_1), path(short_reads_2)

    output:
    tuple val(meta), path("${meta.id}_lr_filtlong.fastq.gz"), emit: reads
    path "versions.yml"                                     , emit: versions

    script:
    """
    filtlong \
        -1 ${short_reads_1} \
        -2 ${short_reads_2} \
        --min_length ${params.longreads_min_length} \
        --keep_percent ${params.longreads_keep_percent} \
        --trim \
        --length_weight ${params.longreads_length_weight} \
        ${long_reads} | gzip > ${meta.id}_lr_filtlong.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filtlong: \$(filtlong --version | sed -e "s/Filtlong v//g")
    END_VERSIONS
    """
}

