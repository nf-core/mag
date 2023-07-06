process NANOPLOT {
    tag "$meta.id"

    conda "bioconda::nanoplot=1.26.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanoplot:1.26.3--py_0' :
        'quay.io/biocontainers/nanoplot:1.26.3--py_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    path '*.png'        , emit: png
    path '*.html'       , emit: html
    path '*.txt'        , emit: txt
    path "versions.yml" , emit: versions

    script:
    def prefix = task.ext.prefix ? "-p ${task.ext.prefix}_" : ''
    def title  = task.ext.prefix ? "${meta.id}_${task.ext.prefix}" : "${meta.id}"
    """
    NanoPlot -t ${task.cpus} \
            ${prefix} \
            --title ${title} \
            -c darkblue \
            --fastq ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        NanoPlot: \$(NanoPlot --version | sed -e "s/NanoPlot //g")
    END_VERSIONS
    """
}
