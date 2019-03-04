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
      --keep_phix                   Keep reads similar to the Illumina internal standard PhiX genome (default: false)

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

    Bin quality check:
      --no_busco                    Disable bin QC with BUSCO (default: false)
      --busco_reference             Download path for BUSCO database, available databases are listed here: https://busco.ezlab.org/
                                    (default: https://busco.ezlab.org/datasets/bacteria_odb9.tar.gz)

    Long read preprocessing:
      --no_nanolyse                 Keep reads similar to the ONT internal standard lambda genome (default: false)
      --filtlong_min_length         Discard any read which is shorter than this value (default: 1000)
      --filtlong_keep_percent       keep the best % of reads (default: 90)
      --filtlong_split              Split reads whenever so much consequence bases fail to match a k-mer in the reference (default: 1000)
      --filtlong_length_weight      The higher the more important is read length when choosing the best reads (default: 10)

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false
params.manifest = false
params.busco_reference = "https://busco.ezlab.org/datasets/bacteria_odb9.tar.gz"

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
params.keep_phix

/*
 * binning options
 */
params.refinem = false
params.no_checkm = false
params.no_busco = false
params.min_contig_size = 1500
params.delta_cont = 5
params.merged_cont = 15
params.delta_compl = 10
params.abs_delta_cov = 0.75
params.delta_gc = 3

params.refinem_db = false
params.ssu_evalue = 1e-6

/*
 * long read preprocessing options
 */
params.no_nanolyse = false
params.filtlong_min_length = 1000
params.filtlong_keep_percent = 90
params.filtlong_split = 1000
params.filtlong_length_weight = 10

/*
 * Create a channel for input read files
 */
if(!params.no_busco){
    Channel
        .fromPath( "${params.busco_reference}", checkIfExists: true )
        .set { file_busco_db }
} else {
    Channel.from()
}

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
        .into { files_sr; files_all_raw }
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
        files_all_raw = Channel.from()
     } else {
         Channel
             .from(params.readPaths)
             .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
             .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
             .into { read_files_fastqc; read_files_fastp }
        files_all_raw = Channel.from()
     }
 } else {
     Channel
         .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
         .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
         .into { read_files_fastqc; read_files_fastp }
    files_all_raw = Channel.from()
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
 * TODO: all new programs!
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
    #chd=\$(readlink -f checkm_data)
    #printf "\$chd\\n\$chd\\n" | checkm data setRoot

    #checkm -h > v_checkm.txt
    echo unknown > v_checkm.txt
    #refinem -h > v_refinem.txt
    echo unknown > v_refinem.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}


/*
 * Trim adapter sequences on long read nanopore files
 */
process porechop { 
    tag "$id"
        
    input:
    set id, lr, sr1, sr2 from files_all_raw
    
    output:
    set id, file("${id}_porechop.fastq"), sr1, sr2 into files_porechop
    set id, lr, val("raw") into files_nanoplot_raw
    
    script:
    """
    porechop -i ${lr} -t "${task.cpus}" -o ${id}_porechop.fastq
    """
}

/*
 * Remove reads mapping to the lambda genome.
 * TODO: add lambda phage to igenomes.config?
 */
nanolyse_reference = "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Escherichia_virus_Lambda/all_assembly_versions/GCA_000840245.1_ViralProj14204/GCA_000840245.1_ViralProj14204_genomic.fna.gz"
if (!params.no_nanolyse) {
    Channel
        .fromPath( "${nanolyse_reference}", checkIfExists: true )
        .set { file_nanolyse_db }
    process nanolyse { 
        tag "$id"

        publishDir "${params.outdir}", mode: 'copy',
            saveAs: {filename -> filename.indexOf(".fastq.gz") == -1 ? "nanolyse/$filename" : null}
            
        input:
        set id, file(lr), sr1, sr2, file(nanolyse_db) from files_porechop.combine(file_nanolyse_db)

        output:
        set id, file("${id}_nanolyse.fastq.gz"), sr1, sr2 into files_nanolyse
        file("${id}_nanolyse_log.txt")
    
        script:
        """
        cat ${lr} | NanoLyse --reference $nanolyse_db | gzip > ${id}_nanolyse.fastq.gz
        
        echo "NanoLyse reference: $nanolyse_reference" >${id}_nanolyse_log.txt
        cat ${lr} | echo "total reads before NanoLyse: \$((`wc -l`/4))" >>${id}_nanolyse_log.txt
        zcat ${id}_nanolyse.fastq.gz | echo "total reads after NanoLyse: \$((`wc -l`/4))" >>${id}_nanolyse_log.txt
        """
    }
} else {
    files_porechop
        .set{ files_nanolyse }
}


/*
 * Quality filter long reads focus on length instead of quality to improve assembly size
 * TODO: Should illumina reads be trimmed/processed before using them here?
 */
process filtlong {
    tag "$id"

    input: 
    set id, lr, sr1, sr2 from files_nanolyse
    
    output:
    set id, file("${id}_lr_filtlong.fastq.gz") into files_lr_filtered 
    set id, file("${id}_lr_filtlong.fastq.gz"), val('filtered') into files_nanoplot_filtered    

    script:
    """
    filtlong \
        -1 ${sr1} \
        -2 ${sr2} \
        --min_length ${params.filtlong_min_length} \
        --keep_percent ${params.filtlong_keep_percent} \
        --trim \
        --split ${params.filtlong_split} \
        --length_weight ${params.filtlong_length_weight} \
        ${lr} | gzip > ${id}_lr_filtlong.fastq.gz
    """
}


/*
 * Quality check for nanopore reads and Quality/Length Plots
 */
process nanoplot {
    tag "$id"
    publishDir "${params.outdir}/${id}/qc/longread_${type}/", mode: 'copy'
    
    input:
    set id, lr, type from files_nanoplot_raw.mix(files_nanoplot_filtered)

    output:
    file '*.png'
    file '*.html'
    file '*.txt'
    
    script:
    """
    NanoPlot -t "${task.cpus}" -p ${type}_  --title ${id}_${type} -c darkblue --fastq ${lr}
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
    set val(name), file("${name}_trimmed*.fastq.gz") into trimmed_reads

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

/*
 * Remove PhiX contamination from Illumina reads
 * TODO: function downloaddb that also makes index
 */
if(!params.keep_phix) {
    phix_reference = "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/Enterobacteria_phage_phiX174_sensu_lato/all_assembly_versions/GCA_002596845.1_ASM259684v1/GCA_002596845.1_ASM259684v1_genomic.fna.gz"
    Channel
        .fromPath( "${phix_reference}", checkIfExists: true )
        .set { file_phix_db }
    process remove_phix {
        tag "$name"

        publishDir "${params.outdir}", mode: 'copy',
            saveAs: {filename -> filename.indexOf(".fastq.gz") == -1 ? "remove_phix/$filename" : null}

        input:
        set val(name), file(reads), file(genome) from trimmed_reads.combine(file_phix_db)

        output:
        set val(name), file("*.fastq.gz") into (trimmed_reads_megahit, trimmed_reads_metabat, trimmed_reads_fastqc, trimmed_sr_spades)
        file("${name}_remove_phix_log.txt")

        script:
        if ( !params.singleEnd ) {
            """
            bowtie2-build --threads "${task.cpus}" "${genome}" ref
            bowtie2 -p "${task.cpus}" -x ref -1 "${reads[0]}" -2 "${reads[1]}" --un-conc-gz unmapped_%.fastq.gz

            echo "Bowtie2 reference: ${genome}" >${name}_remove_phix_log.txt
            zcat ${reads[0]} | echo "Read pairs before removal: \$((`wc -l`/4))" >>${name}_remove_phix_log.txt
            zcat unmapped_1.fastq.gz | echo "Read pairs after removal: \$((`wc -l`/4))" >>${name}_remove_phix_log.txt
            """
        } else {
            """
            bowtie2-build --threads "${task.cpus}" "${genome}" ref
            bowtie2 -p "${task.cpus}" -x ref -U ${reads}  --un-gz unmapped.fastq.gz

            echo "Bowtie2 reference: $ref" >${name}_remove_phix_log.txt
            zcat ${reads[0]} | echo "Reads before removal: \$((`wc -l`/4))" >>${name}_remove_phix_log.txt
            zcat unmapped_1.fastq.gz | echo "Reads after removal: \$((`wc -l`/4))" >>${name}_remove_phix_log.txt
            """
        }

    }
} else {
    trimmed_reads.into {trimmed_reads_megahit; trimmed_reads_metabat; trimmed_reads_fastqc; trimmed_sr_spades}
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
    publishDir "${params.outdir}/", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".fastq.gz") == -1 ? "$filename" : null}

    input:
    set val(name), file(reads) from trimmed_reads_megahit

    output:
    set val(name), file("megahit/${name}.contigs.fa") into (assembly_megahit_to_quast, assembly_megahit_to_refinem)
    set val(name), file("megahit/${name}.contigs.fa"), file(reads) into assembly_megahit_to_metabat

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


/*
 * metaSpades hybrid Assembly
 * TODO: use assembly_graph_spades
 */

 files_lr_filtered
    .combine(trimmed_sr_spades, by: 0)
    .set { files_pre_spades }

process spades {
    tag "$id"
    publishDir "${params.outdir}/${id}/assembly_spades", mode: 'copy', pattern: "${id}*",
        saveAs: {filename -> filename.indexOf(".fastq.gz") == -1 ? filename : null}

    input:
    set id, file(lr), file(sr) from files_pre_spades  

    output:
    //set id, val('spades'), file("${id}_graph_spades.gfa") into assembly_graph_spades
    set val("spades-$id"), file("${id}_scaffolds_spades.fasta") into (assembly_spades_to_quast, assembly_spades_to_refinem)
    set val("spades-$id"), file("${id}_scaffolds_spades.fasta"), file(sr) into assembly_spades_to_metabat
    file("${id}_contigs_spades.fasta")
    file("${id}_spades.log")

    when:
    params.manifest
    !params.singleEnd
     
    script:
    def maxmem = "${task.memory.toString().replaceAll(/[\sGB]/,'')}"
    """
    metaspades.py \
        --threads "${task.cpus}" \
        --memory "$maxmem" \
        --pe1-1 ${sr[0]} \
        --pe1-2 ${sr[1]} \
        --nanopore ${lr} \
        -o spades
    cp spades/assembly_graph_with_scaffolds.gfa ${id}_graph_spades.gfa
    cp spades/scaffolds.fasta ${id}_scaffolds_spades.fasta
    cp spades/contigs.fasta ${id}_contigs_spades.fasta
    cp spades/spades.log ${id}_spades.log
    """
}


process quast {
    tag "$name"
    // publishDir "${params.outdir}/quast/$name", mode: 'copy'

    input:
    set val(name), file(assembly) from assembly_spades_to_quast.mix(assembly_megahit_to_quast)

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
    publishDir "${params.outdir}/metabat", mode: 'copy',
        saveAs: {filename -> 
            if      (filename.indexOf(".fastq.gz") > -1)   filename
            else if (filename.indexOf(".bam") > -1)        filename }

    input:
    set val(name), file(assembly), file(reads) from assembly_spades_to_metabat.mix(assembly_megahit_to_metabat)
    val(min_size) from params.min_contig_size

    output:
    set val(name), file("bins/*") into metabat_bins mode flatten
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

        #if bin foolder is empty
        if [ -z \"\$(ls -A bins)\" ]; then 
            cp ${assembly} bins/
        fi
        """
    } else {
        """
        bowtie2-build --threads "${task.cpus}" "${assembly}" ref
        bowtie2 -p "${task.cpus}" -x ref -U ${reads} | \
            samtools view -@ "${task.cpus}" -bS | \
            samtools sort -@ "${task.cpus}" -o "${name}.bam"
        samtools index "${name}.bam"
        jgi_summarize_bam_contig_depths --outputDepth depth.txt "${name}.bam"
        metabat2 -t "${task.cpus}" -i "${assembly}" -a depth.txt -o bins/"${name}" -m ${min_size}

        #if bin foolder is empty
        if [ -z \"\$(ls -A bins)\" ]; then 
            cp ${assembly} bins/
        fi
        """
    }

}


process busco_download_db {
    tag "${database.baseName}"

    input:
    file(database) from file_busco_db

    output:
    set val("${database.toString().replace(".tar.gz", "")}"), file("buscodb/*") into busco_db

    script:
    """
    mkdir buscodb
    tar -xf ${database} -C buscodb
    """
}

metabat_bins
    .combine(busco_db)
    .set { metabat_db_busco }

/*
 * BUSCO: Quantitative measures for the assessment of genome assembly
 */
process busco {
    tag "${name}-${assembly}"
    publishDir "${params.outdir}/busco", mode: 'copy'

    input:
    set val(name), file(assembly), val(db_name), file(db) from metabat_db_busco

    output:
    file("short_summary_${name}-${assembly}.txt") into busco_summary
    file("${name}-${assembly}_busco_figure.png")
    file("${name}-${assembly}_busco_figure.R")
    file("${name}-${assembly}_busco_log.txt")

    script:
    def binName = "${name}-${assembly}"
    """
    run_BUSCO.py \
        --in ${assembly} \
        --lineage_path $db_name \
        --cpu "${task.cpus}" \
        --blast_single_core \
        --mode genome \
        --out ${binName} \
        >${binName}_busco_log.txt
    generate_plot.py \
        --working_directory run_${binName}
    cp run_${binName}/busco_figure.png ${binName}_busco_figure.png
    cp run_${binName}/busco_figure.R ${binName}_busco_figure.R
    cp run_${binName}/short_summary_${binName}.txt short_summary_${binName}.txt
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
    file ('short_summary_*.txt') from busco_summary.collect()

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
