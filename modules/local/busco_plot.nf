process BUSCO_PLOT {
    tag "${meta.assembler}-${meta.binner}-${meta.id}"

    conda (params.enable_conda ? "bioconda::busco=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.1.0--py_1' :
        'quay.io/biocontainers/busco:5.1.0--py_1' }"

    input:
    tuple val(meta), path(summaries)

    output:
    path("${meta.assembler}-${meta.binner}-${meta.id}.*.busco_figure.png") , optional:true, emit: png
    path("${meta.assembler}-${meta.binner}-${meta.id}.*.busco_figure.R")   , optional:true, emit: rscript
    path "versions.yml"                                                    , emit: versions

    script:
    """
    if [ -n "${summaries}" ]
    then
        # replace dots in bin names within summary file names by underscores
        # currently (BUSCO v5.1.0) generate_plot.py does not allow further dots
        for sum in ${summaries}; do
            if [[ \${sum} =~ short_summary.([_[:alnum:]]+).([_[:alnum:]]+).${meta.assembler}-${meta.binner}-${meta.id}.(.+).txt ]]; then
                mode=\${BASH_REMATCH[1]}
                db_name=\${BASH_REMATCH[2]}
                bin="${meta.assembler}-${meta.binner}-${meta.id}.\${BASH_REMATCH[3]}"
                bin_new="\${bin//./_}"
                mv \${sum} short_summary.\${mode}.\${db_name}.\${bin_new}.txt
            else
                echo "ERROR: the summary filename \${sum} does not match the expected format 'short_summary.([_[:alnum:]]+).([_[:alnum:]]+).${meta.assembler}-${meta.binner}-${meta.id}.(.+).txt'!"
                exit 1
            fi
        done
        generate_plot.py --working_directory .

        mv busco_figure.png "${meta.assembler}-${meta.binner}-${meta.id}.\${mode}.\${db_name}.busco_figure.png"
        mv busco_figure.R "${meta.assembler}-${meta.binner}-${meta.id}.\${mode}.\${db_name}.busco_figure.R"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        R: \$(R --version 2>&1 | sed -n 1p | sed 's/R version //' | sed 's/ (.*//')
        busco: \$(busco --version 2>&1 | sed 's/BUSCO //g')
    END_VERSIONS
    """
}
