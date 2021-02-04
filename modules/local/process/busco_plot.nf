// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process BUSCO_PLOT {
    tag "$assembler-$name"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::busco=4.1.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/busco:4.1.4--py_1"
    } else {
        container "quay.io/biocontainers/busco:4.1.4--py_1"
    }

    input:
    tuple val(assembler), val(name), path(summaries)

    output:
    path("${assembler}-${name}-busco_figure.png")
    path("${assembler}-${name}-busco_figure.R")
    path("${assembler}-${name}-busco_summary.txt")

    script:
    """
    # replace dots in bin names within summary file names by underscores
    # currently (BUSCO v4.1.3) generate_plot.py does not allow further dots
    for sum in ${summaries}; do
        [[ \${sum} =~ short_summary.(.*).${assembler}-${name}.(.*).txt ]];
        db_name=\${BASH_REMATCH[1]}
        bin="${assembler}-${name}.\${BASH_REMATCH[2]}"
        bin_new="\${bin//./_}"
        mv \${sum} short_summary.\${db_name}.\${bin_new}.txt
    done
    generate_plot.py --working_directory .

    mv busco_figure.png ${assembler}-${name}-busco_figure.png
    mv busco_figure.R ${assembler}-${name}-busco_figure.R

    summary_busco.py short_summary.*.txt > ${assembler}-${name}-busco_summary.txt
    """
}
