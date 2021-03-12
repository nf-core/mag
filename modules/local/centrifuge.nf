// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process CENTRIFUGE {
    tag "${meta.id}-${db_name}"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::centrifuge=1.0.4_beta" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/centrifuge:1.0.4_beta--he513fc3_5"
    } else {
        container "quay.io/biocontainers/centrifuge:1.0.4_beta--he513fc3_5"
    }

    input:
    tuple val(meta), path(reads)
    tuple val(db_name), path(db)

    output:
    tuple val("centrifuge"), val(meta), path("results.krona"), emit: results_for_krona
    path "report.txt"
    path "kreport.txt"
    path '*.version.txt'                                     , emit: version

    script:
    def software = getSoftwareName(task.process)
    def input = meta.single_end ? "-U \"${reads}\"" :  "-1 \"${reads[0]}\" -2 \"${reads[1]}\""
    """
    centrifuge -x "${db_name}" \
        -p ${task.cpus} \
        --report-file report.txt \
        -S results.txt \
        $input
    centrifuge-kreport -x "${db_name}" results.txt > kreport.txt
    cat results.txt | cut -f 1,3 > results.krona

    centrifuge --version | head -n 1 | sed 's/^.*centrifuge-class version //' > ${software}.version.txt
    """
}
