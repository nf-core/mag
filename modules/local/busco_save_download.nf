// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process BUSCO_SAVE_DOWNLOAD {

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    input:
    path(busco_downloads)

    output:
    path(busco_downloads, includeInputs: true)

    script:
    """
    """
}
