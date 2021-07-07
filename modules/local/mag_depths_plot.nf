// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process MAG_DEPTHS_PLOT {
    tag "${meta.assembler}-${meta.id}"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:"${meta.assembler}-${meta.id}") }

    conda (params.enable_conda ? "conda-forge::python=3.9 conda-forge::pandas=1.3.0 anaconda::seaborn=0.11.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-d14219255233ee6cacc427e28a7caf8ee42e8c91:0a22c7568e4a509925048454dad9ab37fa8fe776-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-d14219255233ee6cacc427e28a7caf8ee42e8c91:0a22c7568e4a509925048454dad9ab37fa8fe776-0"
    }

    input:
    tuple val(meta), path(depths)
    path(sample_groups)

    output:
    tuple val(meta), path("${meta.assembler}-${meta.id}-binDepths.heatmap.png"), emit: heatmap

    script:
    def software = getSoftwareName(task.process)
    """
    plot_mag_depths.py --bin_depths ${depths} \
                       --groups ${sample_groups} \
                       --out "${meta.assembler}-${meta.id}-binDepths.heatmap.png"
    """
}
