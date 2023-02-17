process TIARA_CLASSIFY {

    conda "conda-forge::r-tidyverse=1.3.1 conda-forge::r-optparse=1.7.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1021c2bc41756fa99bc402f461dad0d1c35358c1:b0c847e4fb89c343b04036e33b2daa19c4152cf5-0' :
        'quay.io/biocontainers/mulled-v2-1021c2bc41756fa99bc402f461dad0d1c35358c1:b0c847e4fb89c343b04036e33b2daa19c4152cf5-0' }"

    input:
    tuple val(meta), path(classification), path(contig2bin)

    output:
    tuple val(meta), path("eukarya/*.fa") , emit: eukarya_bins, optional: true
    tuple val(meta), path("prokarya/*.fa"), emit: prokarya_bins, optional: true

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    domain_classification.R \
        --classification_file ${classification} \
        --contig_to_bin ${contig2bin} \
        ${args} \
        --output_prefix ${prefix}
    """
}
