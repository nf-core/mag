#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/mag
========================================================================================
 nf-core/mag Analysis Pipeline. Started 2018-05-22.
 #### Homepage / Documentation
 https://github.com/nf-core/mag
 #### Authors
 Hadrien Gourlé HadrienG <hadrien.gourle@slu.se> - hadriengourle.com>
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     nf-core/mag v${params.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/mag --reads '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --genome                      Name of iGenomes reference
      -profile                      Hardware config to use. docker / aws

    Options:
      --singleEnd                   Specifies that the input is single end reads

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false

multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

// Validate inputs
if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}
//
// NOTE - THIS IS NOT USED IN THIS PIPELINE, EXAMPLE ONLY
// If you want to use the above in a process, define the following:
//   input:
//   file fasta from fasta
//


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

/*
 * Create channels for input read files
 */
Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { read_files_atropos }


// Header log info
log.info "========================================="
log.info " nf-core/mag v${params.version}"
log.info "========================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Fasta Ref']    = params.fasta
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container']    = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

/*
 * Parse software version numbers
 */
process get_software_versions {
    validExitStatus 0,2

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $params.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    multiqc --version > v_multiqc.txt
    atropos | head -2 | tail -1 > v_atropos.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}


/*
 * STEP 1 - Read trimming and pre/post qc
 */
process atropos {
    tag "$name"
    publishDir "${params.outdir}/atropos", mode: 'copy'

    input:
    set val(name), file(reads) from read_files_atropos

    output:
    set val(name), file("${name}_trimmed_R{1,2}.fastq.gz") into trimmed_reads
    file("${name}.report.txt") into atropos_report

    script:
    def single = reads instanceof Path
    if ( !single ) {
        """
        atropos trim --op-order GACQW --trim-n -o "${name}_trimmed_R1.fastq.gz" \
        -p "${name}_trimmed_R2.fastq.gz" --report-file "${name}.report.txt" \
        --no-cache-adapters --stats both -pe1 "${reads[0]}" -pe2 "${reads[1]}"
        """
    }
    else {
        """
        echo not implemented :-(
        """
    }
}


/*
 * STEP 2 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file ('software_versions/*') from software_versions_yaml.collect()
    file atropos_report from atropos_report

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}
//
//
// /*
//  * STEP 3 - Output Description HTML
//  */
// process output_documentation {
//     tag "$prefix"
//     publishDir "${params.outdir}/Documentation", mode: 'copy'
//
//     input:
//     file output_docs
//
//     output:
//     file "results_description.html"
//
//     script:
//     """
//     markdown_to_html.r $output_docs results_description.html
//     """
// }



/*
 * Completion e-mail notification
 */
 workflow.onComplete {

     // Set up the e-mail variables
     def subject = "[{{ cookiecutter.pipeline_name }}] Successful: $workflow.runName"
     if(!workflow.success){
       subject = "[{{ cookiecutter.pipeline_name }}] FAILED: $workflow.runName"
     }
     def email_fields = [:]
     email_fields['version'] = params.version
     email_fields['runName'] = custom_runName ?: workflow.runName
     email_fields['success'] = workflow.success
     email_fields['dateComplete'] = workflow.complete
     email_fields['duration'] = workflow.duration
     email_fields['exitStatus'] = workflow.exitStatus
     email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
     email_fields['errorReport'] = (workflow.errorReport ?: 'None')
     email_fields['commandLine'] = workflow.commandLine
     email_fields['projectDir'] = workflow.projectDir
     email_fields['summary'] = summary
     email_fields['summary']['Date Started'] = workflow.start
     email_fields['summary']['Date Completed'] = workflow.complete
     email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
     email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
     if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
     if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
     if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
     email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
     email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
     email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

     // Render the TXT template
     def engine = new groovy.text.GStringTemplateEngine()
     def tf = new File("$baseDir/assets/email_template.txt")
     def txt_template = engine.createTemplate(tf).make(email_fields)
     def email_txt = txt_template.toString()

     // Render the HTML template
     def hf = new File("$baseDir/assets/email_template.html")
     def html_template = engine.createTemplate(hf).make(email_fields)
     def email_html = html_template.toString()

     // Render the sendmail template
     def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
     def sf = new File("$baseDir/assets/sendmail_template.txt")
     def sendmail_template = engine.createTemplate(sf).make(smail_fields)
     def sendmail_html = sendmail_template.toString()

     // Send the HTML e-mail
     if (params.email) {
         try {
           if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
           // Try to send HTML e-mail using sendmail
           [ 'sendmail', '-t' ].execute() << sendmail_html
           log.info "[{{ cookiecutter.pipeline_name }}] Sent summary e-mail to $params.email (sendmail)"
         } catch (all) {
           // Catch failures and try with plaintext
           [ 'mail', '-s', subject, params.email ].execute() << email_txt
           log.info "[{{ cookiecutter.pipeline_name }}] Sent summary e-mail to $params.email (mail)"
         }
     }

     // Write summary e-mail HTML to a file
     def output_d = new File( "${params.outdir}/Documentation/" )
     if( !output_d.exists() ) {
       output_d.mkdirs()
     }
     def output_hf = new File( output_d, "pipeline_report.html" )
     output_hf.withWriter { w -> w << email_html }
     def output_tf = new File( output_d, "pipeline_report.txt" )
     output_tf.withWriter { w -> w << email_txt }

     log.info "[{{ cookiecutter.pipeline_name }}] Pipeline Complete"

 }
