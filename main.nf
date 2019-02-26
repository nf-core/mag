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
Daniel Straub <d4straub@gmail.com>
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
    nf-core/mag v${workflow.manifest.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/mag --reads '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Configuration profile to use. Can use multiple (comma separated)
                                    Available: standard, conda, docker, singularity, awsbatch, test

    Hybrid assembly:
      --manifest                    Path to manifest file (must be surrounded with quotes)
                                    Has 4 headerless columns (tab separated): Sample_Id, Long_Reads, Short_Reads_1, Short_Reads_2
                                    Only one file path per entry allowed, join multiple longread files if possible

    Options:
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
      --singleEnd                   Specifies that the input is single end reads
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits

    Trimming options:
      --adapter_forward             Sequence of 3' adapter to remove in the forward reads
      --adapter_reverse             Sequence of 3' adapter to remove in the reverse reads
      --mean_quality                Mean qualified quality value for keeping read (default: 15)
      --trimming_quality            Trimming quality value for the sliding window (default: 15)

    Binning options:
      --refinem                     Enable bin refinement with refinem.
      --refinem_db                  Path to refinem database
      --no_checkm                   Disable bin QC and merging with checkm
      --min_contig_size             Minimum contig size to be considered for binning (default: 1500)
      --delta_cont                  Maximum increase in contamination to merge compatible bins (default: 5)
      --merged_cont                 Maximum total contamination to merge compatible bins (default: 15)
      --delta_compl                 Minimum increase in completion to merge compatible bins (default: 10)
      --abs_delta_cov               Minimum coverage ratio to merge compatible bins (default: 0.75)
      --delta_gc                    Maximum %GC difference to merge compatible bins (default: 3)
      --ssu_evalue                  Evalue threshold to filter incongruent 16S (default 1e-6)

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
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
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false

ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")

// AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}
// Check workDir/outdir paths to be S3 buckets if running on AWSBatch
// related: https://github.com/nextflow-io/nextflow/issues/813
if( workflow.profile == 'awsbatch') {
    if(!workflow.workDir.startsWith('s3:') || !params.outdir.startsWith('s3:')) exit 1, "Workdir or Outdir not on S3 - specify S3 Buckets for each to run on AWSBatch!"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

/*
 * trimming options
 */
params.adapter_forward = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
params.adapter_reverse = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
params.mean_quality = 15
params.trimming_quality = 15

/*
 * binning options
 */
params.refinem = false
params.no_checkm = false
params.min_contig_size = 1500
params.delta_cont = 5
params.merged_cont = 15
params.delta_compl = 10
params.abs_delta_cov = 0.75
params.delta_gc = 3

params.refinem_db = false
params.ssu_evalue = 1e-6

/*
 * Create a channel for input read files
 */
def returnFile(it) {
// Return file if it exists
    if (workflow.profile in ['test', 'localtest'] ) {
        inputFile = file("$workflow.projectDir/" + it)
    } else {
        inputFile = file(it)
    }
    if (!file(inputFile).exists()) exit 1, "Missing file in TSV file: ${inputFile}, see --help for more information"
    return inputFile
}

if(params.manifest){
    manifestFile = file(params.manifest)
    // extracts read files from TSV and distribute into channels
    Channel
        .from(manifestFile)
        .ifEmpty {exit 1, log.info "Cannot find path file ${tsvFile}"}
        .splitCsv(sep:'\t')
        .map { row ->
            def id = row[0]
            def lr = returnFile( row[1] )
            def sr1 = returnFile( row[2] )
            def sr2 = returnFile( row[3] )
            [ id, lr, sr1, sr2 ]
            }
        .into { files_print; files_sr; files_preprocessing }
    // report samples
    files_print
        .subscribe { log.info "\n$it\n" }
    // prepare input for fastqc
    files_sr
        .map { id, lr, sr1, sr2 -> [ id, [ sr1, sr2 ] ] }
        .into { read_files_fastqc; read_files_fastp }
    
} else if(params.readPaths){
     if(params.singleEnd){
         Channel
             .from(params.readPaths)
             .map { row -> [ row[0], [file(row[1][0])]] }
             .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
             .into { read_files_fastqc; read_files_fastp }
     } else {
         Channel
             .from(params.readPaths)
             .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
             .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
             .into { read_files_fastqc; read_files_fastp }
     }
 } else {
     Channel
         .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
         .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
         .into { read_files_fastqc; read_files_fastp }
 }



// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'

nf-core/mag v${workflow.manifest.version}
======================================================="""
def summary = [:]
summary['Pipeline Name']  = 'nf-core/mag'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-mag-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/mag Workflow Summary'
    section_href: 'https://github.com/nf-core/mag'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

/*
 * Parse software version numbers
 */
process get_software_versions {
    validExitStatus 0

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.manifest.nextflowVersion > v_nextflow.txt
    multiqc --version > v_multiqc.txt
    fastqc --version > v_fastqc.txt
    fastp -v 2> v_fastp.txt
    megahit --version > v_megahit.txt
    metabat2 -h 2> v_metabat.txt || true

    # fake checkm data setRoot so we can read checkm version number
    chd=\$(readlink -f checkm_data)
    printf "\$chd\\n\$chd\\n" | checkm data setRoot

    checkm -h > v_checkm.txt
    refinem -h > v_refinem.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}


/*
 * STEP 1 - Read trimming and pre/post qc
 */
process fastqc_raw {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -t "${task.cpus}" -q $reads
    """
}


process fastp {
    tag "$name"
    // publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    set val(name), file(reads) from read_files_fastp
    val adapter from params.adapter_forward
    val adapter_reverse from params.adapter_reverse
    val qual from params.mean_quality
    val trim_qual from params.trimming_quality

    output:
    set val(name), file("${name}_trimmed*.fastq.gz") into (trimmed_reads_megahit, trimmed_reads_metabat, trimmed_reads_fastqc)

    script:
    if ( !params.singleEnd ) {
        """
        fastp -w "${task.cpus}" -q "${qual}" --cut_by_quality5 \
            --cut_by_quality3 --cut_mean_quality "${trim_qual}"\
            --adapter_sequence=${adapter} --adapter_sequence_r2=${adapter_reverse} \
            -i "${reads[0]}" -I "${reads[1]}" \
            -o ${name}_trimmed_R1.fastq.gz -O ${name}_trimmed_R2.fastq.gz
        """
    }
    else {
        """
        fastp -w "${task.cpus}" -q "${qual}" --cut_by_quality5 \
            --cut_by_quality3 --cut_mean_quality "${trim_qual}"\
            --adapter_sequence="${adapter}" --adapter_sequence_r2="${adapter_reverse}" \
            -i ${reads} -o "${name}_trimmed.fastq.gz"
        """
    }
}


process fastqc_trimmed {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    set val(name), file(reads) from trimmed_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results_trimmed

    script:
    """
    fastqc -t "${task.cpus}" -q ${reads}
    """
}


/*
 * STEP 2 - Assembly
 */
process megahit {
    tag "$name"
    publishDir "${params.outdir}/", mode: 'copy'

    input:
    set val(name), file(reads) from trimmed_reads_megahit

    output:
    set val(name), file("megahit/${name}.contigs.fa") into (assembly_quast, assembly_metabat, assembly_refinem)

    script:
    if ( !params.singleEnd ) {
    """
    megahit -t "${task.cpus}" -1 "${reads[0]}" -2 "${reads[1]}" -o megahit \
        --out-prefix "${name}"
    """
    }
    else {
    """
    megahit -t "${task.cpus}" -r ${reads} -o megahit --out-prefix "${name}"
    """
    }
}


process quast {
    tag "$name"
    // publishDir "${params.outdir}/quast/$name", mode: 'copy'

    input:
    set val(name), file(assembly) from assembly_quast

    output:
    file("quast_${name}/*") into quast_results

    script:
    """
    quast.py --threads "${task.cpus}" -l "${name}" "${assembly}" -o "quast_${name}"
    """
}


/*
 * STEP 3 - Binning
 */
process metabat {
    tag "$name"
    publishDir "${params.outdir}/metabat", mode: 'copy'

    input:
    set val(_name), file(assembly) from assembly_metabat
    set val(name), file(reads) from trimmed_reads_metabat
    val(min_size) from params.min_contig_size

    output:
    set val(name), file("bins/") into metabat_bins
    set val(name), file("${name}.bam") into (mapped_reads_checkm, mapped_reads_refinem)

    script:
    if ( !params.singleEnd ) {
    """
    bowtie2-build --threads "${task.cpus}" "${assembly}" ref
    bowtie2 -p "${task.cpus}" -x ref -1 "${reads[0]}" -2 "${reads[1]}" | \
        samtools view -@ "${task.cpus}" -bS | \
        samtools sort -@ "${task.cpus}" -o "${name}.bam"
    samtools index "${name}.bam"
    jgi_summarize_bam_contig_depths --outputDepth depth.txt "${name}.bam"
    metabat2 -t "${task.cpus}" -i "${assembly}" -a depth.txt -o "bins/${name}" -m ${min_size}
    """
    }
    else {
    """
    bowtie2-build --threads "${task.cpus}" "${assembly}" ref
    bowtie2 -p "${task.cpus}" -x ref -U ${reads} | \
        samtools view -@ "${task.cpus}" -bS | \
        samtools sort -@ "${task.cpus}" -o "${name}.bam"
    samtools index "${name}.bam"
    jgi_summarize_bam_contig_depths --outputDepth depth.txt "${name}.bam"
    metabat2 -t "${task.cpus}" -i "${assembly}" -a depth.txt -o bins/"${name}" -m ${min_size}
    """
    }

}


process checkm_download_db {
    output:
        file("checkm_data") into checkm_db

    when:
        params.no_checkm == false

    script:
    """
    mkdir -p checkm_data && \
    cd checkm_data && \
    curl -L -O https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz && \
    tar xzf checkm_data_2015_01_16.tar.gz && \
    cd .. && \
    printf "checkm_data\ncheckm_data\n" | checkm data setRoot
    """
}

process checkm {
    tag "$name"
    publishDir "${params.outdir}/checkm", mode: 'copy'

    input:
    set val(name), file(bins) from metabat_bins
    set val(_name), file(bam) from mapped_reads_checkm
    val(delta_cont) from params.delta_cont
    val(merged_cont) from params.merged_cont
    val(delta_compl) from params.delta_compl
    val(abs_delta_cov) from params.abs_delta_cov
    val(delta_gc) from params.delta_gc
    file("checkm_data") from checkm_db

    output:
    file("${name}_stats/lineage") into checkm_merge_results
    file("${name}_stats/plots") into checkm_merge_plots
    file("${name}_stats/qa.txt") into checkm_merge_qa
    set val(name), file("${name}/") into checkm_merge_bins

    when:
        params.no_checkm == false

    script:
    """
    # re-run setRoot in case checkm has forgotten where the databases are located
    chd=\$(readlink -f checkm_data)
    printf "\$chd\\n\$chd\\n" | checkm data setRoot

    mkdir -p stats
    checkm lineage_wf -t "${task.cpus}" -x fa "${bins}" stats/lineage > stats/qa.txt
    checkm bin_qa_plot -x fa stats/lineage "${bins}" stats/plots

    samtools index "${bam}"
    checkm coverage -t "${task.cpus}" -x fa "${bins}" stats/coverage.txt "${bam}"
    checkm profile stats/coverage.txt > stats/profile.txt
    checkm tree_qa stats/lineage > stats/tree_qa.txt

    samtools view "${bam}" | awk '{if (NR <=1000) print length(\$10)}' >  stats/read_length.txt

    checkm taxon_set domain Bacteria bacteria.ms
    checkm merge -t "${task.cpus}" -x fa --delta_cont "${delta_cont}" --merged_cont "${merged_cont}" \
        bacteria.ms "${bins}" stats/merger/

    mkdir -p ${name}
    mkdir -p ${name}_stats
    merge_bins.py --delta_cont "${delta_cont}" --merged_cont "${merged_cont}" \
        --delta_compl "${delta_compl}" --abs_delta_cov "${abs_delta_cov}" --delta_gc "${delta_gc}" \
        --profile stats/profile.txt --tree stats/tree_qa.txt \
        --length stats/read_length.txt \
        --merger stats/merger/merger.tsv "${bins}" "${name}"
    checkm lineage_wf -t "${task.cpus}" -x fa "${name}" ${name}_stats/lineage > ${name}_stats/qa.txt
    checkm bin_qa_plot -x fa ${name}_stats/lineage "${name}" ${name}_stats/plots
    """
}


process refinem_download_db {
    publishDir "${params.outdir}/db"
    output:
        file("refinem_databases/") into refinem_db

    when:
        params.refinem == true && params.refinem_db ==false && params.no_checkm == false

    script:
    """
    wget https://storage.googleapis.com/mag_refinem_db/refinem_databases.tar.gz
    tar xzf refinem_databases.tar.gz
    """
}

if (params.refinem_db) {
    refinem_db = file(params.refinem_db)
}
mapped_reads_refinem.into {mapped_reads_refinem_1; mapped_reads_refinem_2}
assembly_refinem.into {assembly_refinem_1; assembly_refinem_2}

process refinem_pass_1 {
    tag "${name}"

    input:
    set val(name), file(bins) from checkm_merge_bins
    set val(_name), file(bam) from mapped_reads_refinem_1
    set val(__name), file(assembly) from assembly_refinem_1
    file(refinem_db) from refinem_db

    output:
    set val(name), file("bins_pass_1") into refinem_bins_pass_1

    when:
        params.refinem == true && params.no_checkm == false

    script:
    """
    # filter by GC / tetra
    samtools index "${bam}"
    refinem scaffold_stats -x fa -c "${task.cpus}" "${assembly}" "${bins}" refinem "${bam}"
    refinem outliers --no_plots refinem/scaffold_stats.tsv refinem
    refinem filter_bins -x fa "${bins}" refinem/outliers.tsv new_bins_tmp_1
    # filter by taxonomic assignment
    refinem call_genes -x fa -c "${task.cpus}" new_bins_tmp_1 gene_calls
    refinem taxon_profile -c "${task.cpus}" gene_calls refinem/scaffold_stats.tsv \
        "${refinem_db}/gtdb_r86.dmnd" "${refinem_db}/gtdb_r86_taxonomy.2018-09-25.tsv" \
        taxon_profile
    refinem taxon_filter -c "${task.cpus}" taxon_profile taxon_filter.tsv
    refinem filter_bins -x fa new_bins_tmp_1 taxon_filter.tsv bins_pass_1
    rename 's/.filtered.filtered.fa/.fa/' bins_pass_1/*.filtered.filtered.fa
    """
}

process refinem_pass_2 {
    tag "${name}"
    publishDir "${params.outdir}/", mode: 'copy'

    input:
    set val(name), file(bins) from refinem_bins_pass_1
    set val(_name), file(bam) from mapped_reads_refinem_2
    set val(__name), file(assembly) from assembly_refinem_2
    file(refinem_db) from refinem_db
    val(ssu_evalue) from params.ssu_evalue

    output:
    set val(name), file("refinem") into refinem_bins

    when:
        params.refinem == true && params.no_checkm == false

    script:
    """
    samtools index "${bam}"
    # filter incongruent 16S
    # first we need to re-run taxon profile on the new bin dir
    refinem scaffold_stats -x fa -c "${task.cpus}" "${assembly}" "${bins}" stats "${bam}"
    refinem call_genes -x fa -c "${task.cpus}" "${bins}" gene_calls
    refinem taxon_profile -c "${task.cpus}" gene_calls stats/scaffold_stats.tsv \
        "${refinem_db}/gtdb_r86.dmnd" "${refinem_db}/gtdb_r86_taxonomy.2018-09-25.tsv" \
        taxon_profile
    # then we can identify incongruent ssu
    refinem ssu_erroneous -x fa -c "${task.cpus}" "${bins}" taxon_profile \
        "${refinem_db}/gtdb_r80_ssu" "${refinem_db}/gtdb_r86_taxonomy.2018-09-25.tsv" ssu
    # TODO inspect top-hit and give e-value threshold to filter a contig
    filter_ssu.py --evalue ${ssu_evalue} ssu/ssu_erroneous.tsv ssu/ssu_filtered.tsv 
    refinem filter_bins -x fa "${bins}" ssu/ssu_filtered.tsv refinem
    rename 's/.filtered.fa/.fa/' refinem/*.filtered.fa
    """
}

// TODO next releases:
// good bins channels (from checkm or refinem) + annotation
// multiqc modules for checkm/refinem

/*
 * STEP 4 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config from ch_multiqc_config.collect().ifEmpty([])
    file ('software_versions/*') from software_versions_yaml.collect()
    file (fastqc_raw:'fastqc/*') from fastqc_results.collect().ifEmpty([])
    file (fastqc_trimmed:'fastqc/*') from fastqc_results_trimmed.collect().ifEmpty([])
    file ('quast*/*') from quast_results.collect()

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f ${rtitle} ${rfilename} --config ${multiqc_config} .
    """
}


/*
 * STEP 5 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/mag] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/mag] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
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
          log.info "[nf-core/mag] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/mag] Sent summary e-mail to $params.email (mail)"
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

    log.info "[nf-core/mag] Pipeline Complete"

}
