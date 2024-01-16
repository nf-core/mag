process BUSCO_SAVE_DOWNLOAD {
    // execute sequentially to avoid artefacts when saving files for multiple busco instances
    maxForks 1

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path(busco_downloads)

    output:
    path('busco_downloads/**', includeInputs: true)

    script:
    """
    """
}
