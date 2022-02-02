process BUSCO_SAVE_DOWNLOAD {
    // execute sequentially to avoid artefacts when saving files for multiple busco instances
    maxForks 1

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    path(busco_downloads)

    output:
    path('busco_downloads/**', includeInputs: true)

    script:
    """
    """
}
