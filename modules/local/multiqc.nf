process MULTIQC {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::multiqc=1.12" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.12--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0' }"

    input:
    path multiqc_files
    path mqc_custom_config
    path 'fastqc_raw/*'
    path 'fastqc_trimmed/*'
    path host_removal
    path 'quast*/*'
    path 'bowtie2log/*'
    path short_summary
    path additional

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    script:
    def args = task.ext.args ?: ''
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    read_type = params.single_end ? "--single_end" : ''
    if ( params.host_fasta || params.host_genome ) {
        """
        # get multiqc parsed data for bowtie2
        multiqc -f $custom_config_file *.bowtie2.log
        multiqc_to_custom_tsv.py ${read_type}
        # run multiqc using custom content file instead of original bowtie2 log files
        multiqc -f $custom_config_file --ignore "*.bowtie2.log" .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
        END_VERSIONS
        """
    } else {
        """
        multiqc -f $args .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
        END_VERSIONS
        """
    }
}
