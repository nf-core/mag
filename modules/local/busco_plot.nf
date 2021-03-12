// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process BUSCO_PLOT {
    tag "${meta.assembler}-${meta.id}"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::busco=4.1.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/busco:4.1.4--py_2"
    } else {
        container "quay.io/biocontainers/busco:4.1.4--py_2"
    }

    input:
    tuple val(meta), path(summaries)

    output:
    path("${meta.assembler}-${meta.id}-busco_figure.png")
    path("${meta.assembler}-${meta.id}-busco_figure.R")
    path("${meta.assembler}-${meta.id}-busco_summary.txt")

    script:
    """
    # replace dots in bin names within summary file names by underscores
    # currently (BUSCO v4.1.3) generate_plot.py does not allow further dots
    for sum in ${summaries}; do
        [[ \${sum} =~ short_summary.(.*).${meta.assembler}-${meta.id}.(.*).txt ]];
        db_name=\${BASH_REMATCH[1]}
        bin="${meta.assembler}-${meta.id}.\${BASH_REMATCH[2]}"
        bin_new="\${bin//./_}"
        mv \${sum} short_summary.\${db_name}.\${bin_new}.txt
    done
    generate_plot.py --working_directory .

    mv busco_figure.png "${meta.assembler}-${meta.id}-busco_figure.png"
    mv busco_figure.R "${meta.assembler}-${meta.id}-busco_figure.R"

    summary_busco.py short_summary.*.txt > "${meta.assembler}-${meta.id}-busco_summary.txt"
    """
}
