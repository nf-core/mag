// Subworkflow with functionality specific to the nf-core/mag pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN   } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { samplesheetToList       } from 'plugin/nf-schema'
include { paramsHelp              } from 'plugin/nf-schema'
include { completionEmail         } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary       } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification          } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE   } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    take:
    version // boolean: Display version and exit
    validate_params // boolean: Boolean whether to validate parameters against the schema at runtime
    _monochrome_logs // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir //  string: The output directory where the results will be saved
    input //  string: Path to input samplesheet
    help // boolean: Display help message and exit
    help_full // boolean: Show the full help message
    show_hidden // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1,
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/mag ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/', '')}" }.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/mag/blob/main/CITATIONS.md
"""
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE(
        nextflow_cli_args
    )

    // Note: normally validateInputParameters() goes here, but
    // as we need to use information from samplesheet from the input channel
    // moved it below

    //
    // Create channels from input file provided through params.input and params.assembly_input
    //

    // Validate FASTQ input
    ch_samplesheet = channel.fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map { meta, sr1, sr2, lr ->
            validateInputSamplesheet(meta, sr1, sr2, lr)
        }

    // if coassemble_group or binning_map_mode is set to not 'own', check if all samples in a group have the same platform
    ch_samplesheet
        .map { meta, _sr1, _sr2, _lr -> [meta.group, meta.sr_platform, meta.lr_platform] }
        .groupTuple(by: 0)
        .map { group, sr_platform, lr_platform ->
            def sr_platforms = sr_platform.unique()
            def lr_platforms = lr_platform.unique()
            if (sr_platforms.size() > 1 || lr_platforms.size() > 1) {
                if (params.coassemble_group) {
                    error("[nf-core/mag] ERROR: Multiple short- or long read sequencing platforms found for group ${group}. if --coassemble_group is used, use same platform for all samples in a group.")
                }
                if (params.binning_map_mode != 'own') {
                    error("[nf-core/mag] ERROR: Multiple short- or long read sequencing platforms found for group ${group}. if --binning_map_mode is not 'own', use same platform for all samples in a group.")
                }
            }
        }

    // if binning_map_mode is set to 'all', check if all samples have the same platform
    if (params.binning_map_mode == "all") {
        ch_samplesheet
            .map { meta, _sr1, _sr2, _lr -> meta.sr_platform }
            .unique()
            .toList()
            .map { platforms ->
                if (platforms.size() > 1) {
                    error("[nf-core/mag] ERROR: Multiple short read sequencing platforms (${platforms.join(", ")}) found in samplesheet. Use same platform for all samples when running with binning_map_mode 'all'.")
                }
            }
        ch_samplesheet
            .map { meta, _sr1, _sr2, _lr -> meta.lr_platform }
            .unique()
            .toList()
            .map { platforms ->
                if (platforms.size() > 1) {
                    error("[nf-core/mag] ERROR: Multiple long read sequencing platforms (${platforms.join(", ")}) found in samplesheet. Use same platform for all samples when running with binning_map_mode 'all'.")
                }
            }
    }

    // Prepare FASTQs channel and separate short and long reads and prepare
    ch_raw_short_reads = ch_samplesheet.map { meta, sr1, sr2, _lr ->
        meta.run = meta.run == [] ? "0" : meta.run
        meta.single_end = params.single_end
        if (params.single_end && sr1) {
            return [meta, [sr1]]
        }
        else if (sr1 && sr2) {
            return [meta, [sr1, sr2]]
        }
    }

    ch_raw_long_reads = ch_samplesheet.map { meta, _sr1, _sr2, lr ->
        if (lr) {
            meta.run = meta.run == [] ? "0" : meta.run
            return [meta, lr]
        }
    }

    // Check already if long reads are provided, for later parameter validation
    def hybrid = false
    ch_raw_long_reads
        .map { meta, lr -> [meta.id, lr] }
        .join(ch_raw_long_reads.map { meta, sr1 -> [meta.id, sr1] }, by: 0, remainder: true)
        .map { _id, lr, sr1 ->
            if (lr && sr1) {
                hybrid = true
            }
        }

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters(
        hybrid
    )

    // Validate PRE-ASSEMBLED CONTIG input when supplied
    if (params.assembly_input) {
        ch_input_assemblies = channel.fromList(
            samplesheetToList(params.assembly_input, "${projectDir}/assets/schema_assembly_input.json")
        )
    }

    // Prepare ASSEMBLY input channel
    if (params.assembly_input) {
        ch_input_assemblies.map { meta, fasta ->
            return [meta + [id: params.coassemble_group ? "group-${meta.group}" : meta.id], [fasta]]
        }
    }
    else {
        ch_input_assemblies = channel.empty()
    }

    // Cross validation of input assembly and read IDs: ensure groups are all represented between reads and assemblies
    if (params.assembly_input) {
        ch_read_ids = ch_samplesheet
            .map { meta, _sr1, _sr2, _lr -> params.coassemble_group ? meta.group : meta.id }
            .unique()
            .toList()
            .sort()

        ch_assembly_ids = ch_input_assemblies
            .map { meta, _fasta -> params.coassemble_group ? meta.group : meta.id }
            .unique()
            .toList()
            .sort()

        ch_read_ids
            .concat(ch_assembly_ids)
            .collect(flat: false)
            .map { ids1, ids2 ->
                if (ids1.sort() != ids2.sort()) {
                    exit(1, "[nf-core/mag] ERROR: supplied IDs or Groups in read and assembly CSV files do not match!")
                }
            }
    }

    emit:
    raw_short_reads  = ch_raw_short_reads
    raw_long_reads   = ch_raw_long_reads
    input_assemblies = ch_input_assemblies
    versions         = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {
    take:
    email //  string: email address
    email_on_fail //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url //  string: hook URL for notifications
    multiqc_report //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error("Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters(hybrid) {
    genomeExistsError()

    // Check if binning mapping mode is valid
    if (params.coassemble_group && params.binning_map_mode == 'own') {
        error("[nf-core/mag] ERROR: Invalid combination of parameter '--binning_map_mode own' and parameter '--coassemble_group'. Select either 'all' or 'group' mapping mode when performing group-wise co-assembly.")
    }

    // Check binning length filter parameters are valid
    if (params.bin_max_size && (params.bin_max_size <= params.bin_min_size)) {
        error("[nf-core/mag] ERROR: Invalid value specified for '--bin_max_size'! Value must be greater than --bin_min_size ${params.bin_min_size}. You gave: --bin_max_size ${params.bin_max_size}")
    }

    // Check if settings concerning reproducibility of used tools are consistent and print warning if not
    if (params.megahit_fix_cpu_1 || params.spades_fix_cpus != -1 || params.spadeshybrid_fix_cpus != -1) {
        if (!params.skip_spades && params.spades_fix_cpus == -1) {
            log.warn("[nf-core/mag]: At least one assembly process is run with a parameter to ensure reproducible results, but SPAdes not. Consider using the parameter '--spades_fix_cpus'.")
        }
        if (hybrid && !params.skip_spadeshybrid && params.spadeshybrid_fix_cpus == -1) {
            log.warn("[nf-core/mag]: At least one assembly process is run with a parameter to ensure reproducible results, but SPAdes hybrid not. Consider using the parameter '--spadeshybrid_fix_cpus'.")
        }
        if (!params.skip_megahit && !params.megahit_fix_cpu_1) {
            log.warn("[nf-core/mag]: At least one assembly process is run with a parameter to ensure reproducible results, but MEGAHIT not. Consider using the parameter '--megahit_fix_cpu_1'.")
        }
        if (!params.skip_binning && params.metabat_rng_seed == 0) {
            log.warn("[nf-core/mag]: At least one assembly process is run with a parameter to ensure reproducible results, but for MetaBAT2 a random seed is specified ('--metabat_rng_seed 0'). Consider specifying a positive seed instead.")
        }
    }

    // Check if SPAdes and single_end
    if ((!params.skip_spades || !params.skip_spadeshybrid) && params.single_end) {
        log.warn('[nf-core/mag]: metaSPAdes does not support single-end data. SPAdes will be skipped.')
    }

    // Check if parameters for host contamination removal are valid
    if (params.host_fasta && params.host_genome) {
        error('[nf-core/mag] ERROR: Both host fasta reference and iGenomes genome are specified to remove host contamination! Invalid combination, please specify either --host_fasta or --host_genome.')
    }
    if (hybrid && (params.host_fasta || params.host_genome) && params.longread_filtering_tool == "filtlong" && params.longreads_length_weight > 0) {
        log.warn("[nf-core/mag]: The parameter --longreads_length_weight is ${params.longreads_length_weight}, causing the read length being more important for long read filtering than the read quality. Set --longreads_length_weight to 1 in order to assign equal weights.")
    }
    if (params.host_genome) {
        if (!params.genomes) {
            error('[nf-core/mag] ERROR: No config file containing genomes provided!')
        }
        // Check if host genome exists in the config file
        if (!params.genomes.containsKey(params.host_genome)) {
            error(
                '=============================================================================\n' + "  Host genome '${params.host_genome}' not found in any config files provided to the pipeline.\n" + '  Currently, the available genome keys are:\n' + "  ${params.genomes.keySet().join(', ')}\n" + '==================================================================================='
            )
        }
        if (!params.genomes[params.host_genome].fasta) {
            error("[nf-core/mag] ERROR: No fasta file specified for the host genome ${params.host_genome}!")
        }
        if (!params.genomes[params.host_genome].bowtie2) {
            error("[nf-core/mag] ERROR: No Bowtie 2 index file specified for the host genome ${params.host_genome}!")
        }
    }

    // Check MetaBAT2 inputs
    if (!params.skip_metabat2 && params.min_contig_size < 1500) {
        log.warn("[nf-core/mag]: Specified min. contig size under minimum for MetaBAT2. MetaBAT2 will be run with 1500 (other binners not affected). You supplied: --min_contig_size ${params.min_contig_size}")
    }

    // Check more than one binner is run for bin refinement  (required DAS by Tool)
    // If the number of run binners (i.e., number of not-skipped) is more than one, otherwise throw an error
    if (params.refine_bins_dastool && !([params.skip_metabat2, params.skip_maxbin2, params.skip_concoct, params.skip_metabinner, params.skip_comebin, params.skip_semibin].count(false) > 1)) {
        error('[nf-core/mag] ERROR: Bin refinement with --refine_bins_dastool requires at least two binners to be running (not skipped). Check input.')
    }

    // Check that bin refinement is actually turned on if any of the refined bins are requested for downstream
    if (!params.refine_bins_dastool && params.postbinning_input != 'raw_bins_only') {
        error("[nf-core/mag] ERROR: The parameter '--postbinning_input ${params.postbinning_input}' for downstream steps can only be specified if bin refinement is activated with --refine_bins_dastool! Check input.")
    }

    if (params.skip_binqc && (params.run_busco || params.run_checkm || params.run_checkm2)) {
        error("[nf-core/mag] ERROR: Both --skip_binqc and --run_<bin_qc_tool_name> are specified! Invalid combination, please specify either --skip_binqc or --run_<bin_qc_tool_name>.")
    }
    if (!params.skip_binqc && !params.run_busco && !params.run_checkm && !params.run_checkm2) {
        log.warn('[nf-core/mag]: --skip_binqc is not specified, but no bin quality assessment tool is set to run! Bin QC will be skipped.')
    }

    if (!params.skip_binqc && params.run_busco) {
        if (params.busco_db && !params.busco_db_lineage) {
            log.warn('[nf-core/mag]: WARNING: You have supplied a database to --busco_db - BUSCO will run in offline mode. Please note that BUSCO may fail if you have an incomplete database and are running with --busco_db_lineage auto!')
        }

        if (params.busco_db && file(params.busco_db).isDirectory() && !file(params.busco_db).listFiles().any { file -> file.toString().contains('lineages') }) {
            error("[nf-core/mag] ERROR: Directory supplied to `--busco_db` must contain a `lineages/` subdirectory that itself contains one or more BUSCO lineage files! Check: --busco_db ${params.busco_db}")
        }
    }

    if (!params.skip_gtdbtk) {
        if (params.skip_binqc) {
            log.warn('[nf-core/mag]: --skip_binqc is specified, but --skip_gtdbtk is explicitly set to run! GTDB-tk will be omitted because GTDB-tk bin classification requires bin filtering based on BUSCO or CheckM QC results to avoid GTDB-tk errors.')
        }

        if (!params.run_busco && !params.run_checkm && !params.run_checkm2) {
            log.warn('[nf-core/mag]: GTDB-tk requires bin quality information from BUSCO, CheckM or CheckM2. Please enable at least one of these tools, otherwise GTDB-tk will be skipped.')
        }
    }

    // Check if CAT parameters are valid
    if (params.cat_db && params.cat_db_generate) {
        error('[nf-core/mag] ERROR: Invalid combination of parameters --cat_db and --cat_db_generate is specified! Please specify either --cat_db or --cat_db_generate.')
    }
    if (params.save_cat_db && !params.cat_db_generate) {
        error('[nf-core/mag] ERROR: Invalid parameter combination: parameter --save_cat_db specified, but not --cat_db_generate! Note also that the parameter --save_cat_db does not work in combination with --cat_db.')
    }

    // Check MetaEuk db parameters
    if (params.metaeuk_mmseqs_db && params.metaeuk_db) {
        error('[nf-core/mag] ERROR: Invalid parameter combination: both --metaeuk_mmseqs_db and --metaeuk_db are specified! Please specify either --metaeuk_mmseqs_db or --metaeuk_db.')
    }
    if (params.save_mmseqs_db && !params.metaeuk_mmseqs_db) {
        error('[nf-core/mag] ERROR: Invalid parameter combination: --save_mmseqs_db supplied but no database has been requested for download with --metaeuk_mmseqs_db!')
    }

    // Check Prokka parameters
    if (params.prokka_with_compliance && !params.prokka_compliance_centre) {
        error('[nf-core/mag] ERROR: Invalid parameter combination: running PROKKA with compliance mode requires a centre name specified with `--prokka_compliance_centre <XYZ>`!')
    }

    // Check BIgMAG parameters
    if (params.generate_bigmag_file && (!params.run_gunc || !params.run_checkm2 || !params.run_busco || params.skip_gtdbtk || params.skip_quast || params.skip_binqc)) {
        error('[nf-core/mag] ERROR: To generate the BIgMAG file you need to include the parameters `--run_checkm2` and `--run_gunc`, and you cannot skip BINQC, GTDB-TK, QUAST nor BUSCO.')
    }

    // Check ancient DNA damage parameters
    if (params.ancient_dna && params.binning_map_mode != 'own') {
        log.warn("[nf-core/mag] WARNING: Running in --binning_map_mode ${params.binning_map_mode} will result in unstable pyDamage output files. You might not receive pyDamage results for all bins in bin_summary.tsv, and `-resume` may not work; `--binning_map_mode own` is recommended!")
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(meta, sr1, sr2, lr) {

    if ((!sr2 && !lr) && !params.single_end) {
        error("[nf-core/mag] ERROR: Single-end data must be executed with `--single_end`. Note that it is not possible to mix single- and paired-end data in one run! Check input CSV for sample: ${meta.id}")
    }
    if (sr2 && params.single_end) {
        error("[nf-core/mag] ERROR: Paired-end data must be executed without `--single_end`. Note that it is not possible to mix single- and paired-end data in one run! Check input CSV for sample: ${meta.id}")
    }

    return [meta, sr1, sr2, lr]
}

//
// Get attribute from genome config file e.g. fasta
//
// Note: user uses --host_genome in mag
def getGenomeAttribute(attribute) {
    if (params.genomes && params.host_genome && params.genomes.containsKey(params.host_genome)) {
        if (params.genomes[params.host_genome].containsKey(attribute)) {
            return params.genomes[params.host_genome][attribute]
        }
    }
    return null
}
//
// Exit pipeline if incorrect --genome key provided
//
// Note: user uses --host_genome in mag

def genomeExistsError() {
    if (params.genomes && params.host_genome && !params.genomes.containsKey(params.host_genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" + "  Genome '${params.host_genome}' not found in any config files provided to the pipeline.\n" + "  Currently, the available genome keys are:\n" + "  ${params.genomes.keySet().join(", ")}\n" + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def text_seq_qc = "Sequencing quality control was performed with FastQC (Andrews 2010)."

    def shortread_qc_tools = [
        params.clip_tool == 'fastp' ? "fastp (Chen et al. 2018)" : "",
        params.clip_tool == 'adapterremoval' ? "AdapterRemoval (Schubert et al. 2016)" : "",
        params.clip_tool == 'trimmomatic' ? "Trimmomatic (Bolger et al. 2014)" : "",
    ].findAll { tool -> tool != '' }
    def text_shortread_qc = "Short read preprocessing was performed with ${shortread_qc_tools.join(', ')}."

    def text_mapping = "Read alignment was performed with Bowtie2 (Langmead and Salzberg 2012) and minimap2 (Li 2018)."

    def longread_qc_tools = [
        !params.skip_adapter_trimming && params.longread_adaptertrimming_tool == 'porechop' ? "Porechop (Wick et al. 2017)" : "",
        !params.skip_adapter_trimming && params.longread_adaptertrimming_tool == 'porechop_abi' ? "Porechop ABI (Bonenfant et al. 2022)" : "",
        !params.skip_longread_filtering && params.longread_filtering_tool == 'filtlong' ? "Filtlong (Wick 2019)" : "",
        !params.skip_longread_filtering && params.longread_filtering_tool == 'nanoq' ? "Nanoq (Steinig et al. 2022)" : "",
        !params.skip_longread_filtering && params.longread_filtering_tool == 'chopper' ? "Chopper (De Coster et al. 2018)" : "",
    ].findAll { tool -> tool != '' }
    def text_longread_qc = "Long read preprocessing was performed with ${longread_qc_tools.join(', ')}."

    def text_bbnorm = "Read depth normalisation was carried out with BBNorm (Bushnell et al.) and Seqtk (Li et al.)."

    def assembly_tools = [
        !params.skip_megahit ? "MEGAHIT (Li et al. 2016)" : "",
        !params.skip_spades ? "metaSPAdes (Nurk et al. 2017)" : "",
        !params.skip_spadeshybrid ? "hybridSPAdes (Antipov et al. 2016)" : "",
        !params.skip_metamdbg ? "metaMDBG (Benoit et al. 2024)" : "",
        !params.skip_flye ? "metaFlye (Kolmogorov et al. 2020)" : "",
    ].findAll { tool -> tool != '' }
    def text_assembly = "Metagenome assembly was performed with ${assembly_tools.join(', ')}."

    def assembly_qc_tools = [
        !params.skip_quast ? "metaQUAST (Mikheenko et al. 2016)" : "",
        !params.skip_ale ? "ALE (Clark et al. 2013)" : "",
    ].findAll { tool -> tool != '' }
    def text_assembly_qc = "Assembly quality was assessed with ${assembly_qc_tools.join(', ')}."

    def gene_prediction_tools = [
        !params.skip_prodigal ? "Prodigal (Hyatt et al. 2010)" : "",
        !params.skip_prokka ? "Prokka (Seemann 2014)" : "",
        !params.skip_metaeuk ? "MetaEuk (Levy Karin et al. 2020) with MMseqs2 (Steinegger and Söding 2017)" : "",
    ].findAll { tool -> tool != '' }
    def text_gene_prediction = "Gene prediction was performed with ${gene_prediction_tools.join(', ')}."

    def text_virus_id = "Viral sequence identification was carried out with geNomad (Camargo et al. 2023)."

    def text_tiara = "Eukaryotic contig classification was performed with Tiara (Karlicki et al. 2022)."

    def binning_tools = [
        !params.skip_metabat2 ? "MetaBAT2 (Kang et al. 2019)" : "",
        !params.skip_maxbin2 ? "MaxBin2 (Wu et al. 2015)" : "",
        !params.skip_concoct ? "CONCOCT (Alneberg et al. 2014)" : "",
        !params.skip_comebin ? "COMEBin (Wang et al. 2024)" : "",
        !params.skip_metabinner ? "MetaBinner (Wang et al. 2023)" : "",
        !params.skip_semibin ? "SemiBin2 (Pan et al. 2022)" : "",
    ].findAll { tool -> tool != '' }
    def text_binning = "Metagenome binning was performed with ${binning_tools.join(', ')}."

    def text_bin_refinement = "Bin refinement was performed with DAS Tool (Sieber et al. 2018)."

    def binqc_tools = [
        params.run_busco ? "BUSCO (Seppey et al. 2019)" : "",
        params.run_checkm ? "CheckM (Parks et al. 2015)" : "",
        params.run_checkm2 ? "CheckM2 (Chklovski et al. 2023)" : "",
        params.run_gunc ? "GUNC (Orakov et al. 2021)" : "",
    ].findAll { tool -> tool != '' }
    def text_binqc = "Bin quality assessment was carried out with ${binqc_tools.join(', ')}."

    def text_gtdbtk = "Taxonomic classification of bins was performed with GTDB-Tk (Chaumeil et al. 2020)."

    def text_ancient_dna = [
        "Ancient DNA damage assessment was performed with PyDamage (Borry et al. 2021).",
        !params.skip_ancient_damagecorrection ? "Damage correction was carried out with variant calling using FreeBayes (Garrison and Marth 2012), BCFtools (Danecek et al. 2021), and SAMtools (Li et al. 2009)." : "",
    ]

    def citation_text = [
        "Tools used in the workflow included:",
        text_seq_qc,
        (!params.skip_shortread_qc && !params.skip_clipping) ? text_shortread_qc : "",
        (params.host_fasta || params.host_genome || !params.skip_binning || params.ancient_dna || !params.skip_ale) ? text_mapping : "",
        (!params.skip_longread_qc && longread_qc_tools) ? text_longread_qc : "",
        params.bbnorm ? text_bbnorm : "",
        assembly_tools ? text_assembly : "",
        assembly_qc_tools ? text_assembly_qc : "",
        gene_prediction_tools ? text_gene_prediction : "",
        params.run_virus_identification ? text_virus_id : "",
        !params.skip_tiara ? text_tiara : "",
        (!params.skip_binning && binning_tools) ? text_binning : "",
        (!params.skip_binqc && params.refine_bins_dastool) ? text_bin_refinement : "",
        (!params.skip_binqc && binqc_tools) ? text_binqc : "",
        !params.skip_gtdbtk ? text_gtdbtk : "",
        params.ancient_dna ? text_ancient_dna : "",
        "Pipeline results statistics were summarised with MultiQC (Ewels et al. 2016).",
    ].flatten().join(' ').trim().replaceAll("\\s+", " ")

    return citation_text
}

def toolBibliographyText() {
    def references = [
        "<li>Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data [Online]. URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</li>",
        "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047-3048. doi: 10.1093/bioinformatics/btw354</li>",
    ]

    if (!params.skip_shortread_qc && !params.skip_clipping) {
        if (params.clip_tool == 'fastp') {
            references << "<li>Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884–i890. doi: 10.1093/bioinformatics/bty560</li>"
        }
        else if (params.clip_tool == 'adapterremoval') {
            references << "<li>Schubert, M., Lindgreen, S., & Orlando, L. (2016). AdapterRemoval v2: rapid adapter trimming, identification, and read merging. BMC Research Notes, 9, 88. doi: 10.1186/s13104-016-1900-2</li>"
        }
        else if (params.clip_tool == 'trimmomatic') {
            references << "<li>Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114-2120. doi: 10.1093/bioinformatics/btu170</li>"
        }
    }
    if (params.host_fasta || params.host_genome || !params.skip_binning || params.ancient_dna || !params.skip_ale) {
        references << "<li>Langmead, B. and Salzberg, S. L. 2012 Fast gapped-read alignment with Bowtie 2. Nature methods, 9(4), p. 357–359. doi: 10.1038/nmeth.1923.</li>"
        // Note: we don't have a simple way to determine if long reads are present, so we add minimap2 at the same time as Bowtie2
        references << "<li>Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics , 34(18), 3094–3100. doi: 10.1093/bioinformatics/bty191</li>"
    }
    if (!params.skip_longread_qc && !params.skip_adapter_trimming) {
        if (params.longread_adaptertrimming_tool == 'porechop') {
            references << "<li>Wick RR. (2017). Porechop. URL: https://github.com/rrwick/Porechop</li>"
        }
        else if (params.longread_adaptertrimming_tool == 'porechop_abi') {
            references << "<li>Bonenfant, Q., Noé, L., & Touzet, H. (2022). Porechop_ABI: discovering unknown adapters in ONT sequencing reads for downstream trimming. bioRxiv. 10.1101/2022.07.07.499093</li>"
        }
    }
    if (!params.skip_longread_qc && !params.skip_longread_filtering) {
        if (params.longread_filtering_tool == 'filtlong') {
            references << "<li>Wick RR. (2019). Filtlong. URL: https://github.com/rrwick/Filtlong</li>"
        }
        else if (params.longread_filtering_tool == 'nanoq') {
            references << "<li>Steinig, E., Coin, L. (2022). Nanoq: ultra-fast quality control for nanopore reads. Journal of Open Source Software, 7(69), 2991, doi: 10.21105/joss.02991</li>"
        }
        else if (params.longread_filtering_tool == 'chopper') {
            references << "<li>De Coster W, D'Hert S, Schultz DT, Cruts M, Van Broeckhoven C. (2018). NanoPack: visualizing and processing long-read sequencing data. Bioinformatics. 2018 Aug 1;34(15):2666-2669. doi: 10.1093/bioinformatics/bty149</li>"
        }
    }
    if (params.bbnorm) {
        references << "<li>BBnorm/BBTools. URL: http://sourceforge.net/projects/bbmap/</li>"
        references << "<li>Seqtk. URL: https://github.com/lh3/seqtk</li>"
    }
    if (!params.skip_megahit) {
        references << "<li>Li, D., Luo, R., Liu, C. M., Leung, C. M., Ting, H. F., Sadakane, K., ... & Lam, T. W. (2016). MEGAHIT v1. 0: a fast and scalable metagenome assembler driven by advanced methodologies and community practices. Methods, 102, 3-11. doi: 10.1016/j.ymeth.2016.02.020.</li>"
    }
    if (!params.skip_spades) {
        references << "<li>Nurk, S., Meleshko, D., Korobeynikov, A., & Pevzner, P. A. (2017). metaSPAdes: a new versatile metagenomic assembler. Genome research, 27(5), 824-834. doi: 10.1101/gr.213959.116.</li>"
    }
    if (!params.skip_spadeshybrid) {
        references << "<li>Antipov, D., Korobeynikov, A., McLean, J. S., & Pevzner, P. A. (2016). hybridSPAdes: an algorithm for hybrid assembly of short and long reads. Bioinformatics, 32(7), 1009-1015. doi: 10.1093/bioinformatics/btv688</li>"
    }
    if (!params.skip_metamdbg) {
        references << "<li>Benoit, G., Raguideau, S., James, R. et al. (2024). High-quality metagenome assembly from long accurate reads with metaMDBG. Nat Biotechnol 42, 1378–1383. doi: 10.1038/s41587-023-01983-6</li>"
    }
    if (!params.skip_flye) {
        references << "<li>Kolmogorov, M., Bickhart, D. M., Behsaz, B. et al. (2020). metaFlye: scalable long-read metagenome assembly using repeat graphs. Nature Methods, 17(11), 1103-1110. doi: 10.1038/s41592-020-00971-x</li>"
    }
    if (!params.skip_quast) {
        references << "<li>Mikheenko, A., Saveliev, V., & Gurevich, A. (2016). MetaQUAST: evaluation of metagenome assemblies. Bioinformatics, 32(7), 1088-1090. doi: 10.1093/bioinformatics/btv697</li>"
    }
    if (!params.skip_ale) {
        references << "<li>Clark, S. C., Egan, R., Frazier, P. I., & Wang, Z. (2013). ALE: a generic assembly likelihood evaluation framework for assessing the accuracy of genome and metagenome assemblies. Bioinformatics, 29(4), 435-443. doi: 10.1093/bioinformatics/bts723</li>"
    }
    if (!params.skip_prodigal) {
        references << "<li>Hyatt, D., Chen, G. L., Locascio, P. F., Land, M. L., Larimer, F. W., & Hauser, L. J. (2010). Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics, 11, 119. doi: 10.1186/1471-2105-11-119</li>"
    }
    if (!params.skip_prokka) {
        references << "<li>Seemann, T. (2014). Prokka: rapid prokaryotic genome annotation. Bioinformatics, 30(14), 2068-2069. doi: 10.1093/bioinformatics/btu153</li>"
    }
    if (!params.skip_metaeuk) {
        references << "<li>Levy Karin, E., Mirdita, M. & Söding, J. (2020). MetaEuk—sensitive, high-throughput gene discovery, and annotation for large-scale eukaryotic metagenomics. Microbiome 8, 48. doi: 10.1186/s40168-020-00808-x</li>"
        references << "<li>Steinegger, M., & Söding, J. (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat Biotechnol 35, 1026–1028. doi: 10.1038/nbt.3988</li>"
    }
    if (params.run_virus_identification) {
        references << "<li>Camargo, A. P., et al. (2023). Identification of mobile genetic elements with geNomad. Nature Biotechnology 42, 1303–1312. doi: 10.1038/s41587-023-01953-y</li>"
    }
    if (!params.skip_tiara) {
        references << "<li>Karlicki, M., Antonowicz, S., & Karnkowska, A. (2022). Tiara: deep learning-based classification system for eukaryotic sequences. Bioinformatics, 38(2), 344–350. doi: 10.1093/bioinformatics/btab672</li>"
    }
    if (!params.skip_binning) {
        if (!params.skip_metabat2) {
            references << "<li>Kang, D. D., Li, F., Kirton, E., Thomas, A., Egan, R., An, H., & Wang, Z. (2019). MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. PeerJ, 7, e7359. doi: 10.7717/peerj.7359</li>"
        }
        if (!params.skip_maxbin2) {
            references << "<li>Wu, Y. W., Simmons, B. A., & Singer, S. W. (2015). MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. Bioinformatics, 32(4), 605-607. doi: 10.1093/bioinformatics/btv638</li>"
        }
        if (!params.skip_concoct) {
            references << "<li>Alneberg, J., Bjarnason, B. S., de Bruijn, I., Schirmer, M., Quick, J., Ijaz, U. Z., Lahti, L., Loman, N. J., Andersson, A. F., & Quince, C. (2014). Binning metagenomic contigs by coverage and composition. Nature Methods, 11(11), 1144–1146. doi: 10.1038/nmeth.3103</li>"
        }
        if (!params.skip_comebin) {
            references << "<li>Wang, Z., You, R., Han, H., Liu, W., Sun, F., & Zhu, S. (2024). Effective binning of metagenomic contigs using contrastive multi-view representation learning. Nat Commun. 2024 Jan 17;15(1):585. doi: 10.1038/s41467-023-44290-z</li>"
        }
        if (!params.skip_metabinner) {
            references << "<li>Wang, Z., Huang, P., You, R., Sun, F., & Zhu, S. (2023). MetaBinner: a high-performance and stand-alone ensemble binning method to recover individual genomes from complex microbial communities. Genome Biol. 2023 Jan 6;24(1):1. doi: 10.1186/s13059-022-02832-6</li>"
        }
        if (!params.skip_semibin) {
            references << "<li>Pan, S., Zhu, C., Zhao, X. M., & Coelho, L. P. (2022). A deep siamese neural network improves metagenome-assembled genomes in microbiome datasets across different environments. Nat Commun. 2022 Apr 28;13(1):2326. doi: 10.1038/s41467-022-29843-y</li>"
        }
    }
    if (!params.skip_binqc) {
        if (params.run_busco) {
            references << "<li>Seppey, M., Manni, M., & Zdobnov, E. M. (2019). BUSCO: assessing genome assembly and annotation completeness. In Gene prediction (pp. 227-245). Humana, New York, NY. doi: 10.1007/978-1-4939-9173-0_14</li>"
        }
        if (params.run_checkm) {
            references << "<li>Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P., & Tyson, G. W. (2015). CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Research, 25(7), 1043–1055. doi: 10.1101/gr.186072.114</li>"
        }
        if (params.run_checkm2) {
            references << "<li>Chklovski, A., Parks, D. H., Woodcroft, B. J., & Tyson, G. W. (2023). CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. Nature Methods, 20(8), 1203-1212. doi: 10.1038/s41592-023-01940-w</li>"
        }
        if (params.refine_bins_dastool) {
            references << "<li>Sieber, C. M. K., et al. (2018). Recovery of Genomes from Metagenomes via a Dereplication, Aggregation and Scoring Strategy. Nature Microbiology 3 (7): 836-43. doi: 10.1038/s41564-018-0171-1</li>"
        }
        if (params.run_gunc) {
            references << "<li>Orakov, A., Fullam, A., Coelho, A. P., Khedkar, S., Szklarczyk, D., Mende, D. R., Schmidt, T. S. B., & Bork, P. (2021). GUNC: detection of chimerism and contamination in prokaryotic genomes. Genome Biology, 22(1), 178. doi: 10.1186/s13059-021-02393-0</li>"
        }
    }
    if (!params.skip_gtdbtk) {
        references << "<li>Chaumeil, P. A., Mussig, A. J., Hugenholtz, P., & Parks, D. H. (2020). GTDB-Tk: a toolkit to classify genomes with the Genome Taxonomy Database. Bioinformatics, 36(6), 1925–1927. doi: 10.1093/bioinformatics/btz848</li>"
    }
    if (params.ancient_dna) {
        references << "<li>Borry, M., Hübner, A., Rohrlach, A. B., & Warinner, C. (2021). PyDamage: automated ancient damage identification and estimation for contigs in ancient DNA de novo assembly. PeerJ, 9, e11845. doi: 10.7717/peerj.11845</li>"
        if (!params.skip_ancient_damagecorrection) {
            references << "<li>Garrison, E., & Marth, G. (2012). Haplotype-based variant detection from short-read sequencing. arXiv preprint arXiv:1207.3907 [q-bio.GN]</li>"
            references << "<li>Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Stackmayer, V., Li, H., ... & SAMtools Subgroup. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. doi: 10.1093/gigascience/giab008</li>"
            references << "<li>Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., … 1000 Genome Project Data Processing Subgroup. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. doi: 10.1093/bioinformatics/btp352</li>"
        }
    }

    def reference_text = references.join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    }
    else {
        meta["doi_text"] = ""
    }
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
