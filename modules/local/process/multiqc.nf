// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MULTIQC {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename: filename, options: params.options, publish_dir: getSoftwareName(task.process), publish_id: '') }

    conda (params.enable_conda ? "bioconda::multiqc=1.9" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/multiqc:1.9--pyh9f0ad1d_0"
    } else {
        container "quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    }

    input:
    path multiqc_files
    path mqc_custom_config
    path 'fastqc_raw/*'
    path 'fastqc_trimmed/*'
    path host_removal
    path 'quast*/*'
    path 'bowtie2log/*'
    path short_summary

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "*.version.txt"       , emit: version

    script:
    def software = getSoftwareName(task.process)
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    read_type = params.single_end ? "--single_end" : ''
    if ( params.host_fasta || params.host_genome ) {
        """
        # get multiqc parsed data for bowtie2
        multiqc -f $custom_config_file *.bowtie2.log
        multiqc_to_custom_tsv.py ${read_type}
        # run multiqc using custom content file instead of original bowtie2 log files
        multiqc -f $custom_config_file --ignore "*.bowtie2.log" .
        multiqc --version | sed -e "s/multiqc, version //g" > ${software}.version.txt
        """
    } else {
        """
        multiqc -f $options.args .
        multiqc --version | sed -e "s/multiqc, version //g" > ${software}.version.txt
        """
    }
}
