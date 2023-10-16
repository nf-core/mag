process CENTRIFUGE {
    tag "${meta.id}-${db_name}"

    conda "bioconda::centrifuge=1.0.4_beta"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/centrifuge:1.0.4_beta--he513fc3_5' :
        'biocontainers/centrifuge:1.0.4_beta--he513fc3_5' }"

    input:
    tuple val(meta), path(reads)
    tuple val(db_name), path(db)

    output:
    tuple val("centrifuge"), val(meta), path("results.krona"), emit: results_for_krona
    path "report.txt"                                        , emit: report
    tuple val(meta), path("*kreport.txt")                    , emit: kreport
    path "versions.yml"                                      , emit: versions

    script:
    def input = meta.single_end ? "-U \"${reads}\"" :  "-1 \"${reads[0]}\" -2 \"${reads[1]}\""
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    centrifuge -x "${db_name}" \
        -p ${task.cpus} \
        --report-file report.txt \
        -S results.txt \
        $input
    centrifuge-kreport -x "${db_name}" results.txt > ${prefix}.centrifuge_kreport.txt
    cat results.txt | cut -f 1,3 > results.krona

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        centrifuge: \$(centrifuge --version | sed -n 1p | sed 's/^.*centrifuge-class version //')
    END_VERSIONS
    """
}
