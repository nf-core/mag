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
    def tools = [
        "FastQC (Andrews 2010)",
        "MultiQC (Ewels et al. 2016)",
    ]

    if (!params.skip_shortread_qc && !params.skip_clipping) {
        if (params.clip_tool == 'fastp') {
            tools << "fastp (Chen et al. 2018)"
        }
        else if (params.clip_tool == 'adapterremoval') {
            tools << "AdapterRemoval (Schubert et al. 2016)"
        }
        else if (params.clip_tool == 'trimmomatic') {
            tools << "Trimmomatic (Bolger et al. 2014)"
        }
    }
    if (!params.skip_longread_qc && !params.skip_adapter_trimming) {
        if (params.longread_adaptertrimming_tool == 'porechop') {
            tools << "Porechop (Wick et al. 2017)"
        }
        else if (params.longread_adaptertrimming_tool == 'porechop_abi') {
            tools << "Porechop ABI (Wick et al. 2017)"
        }
    }
    if (!params.skip_longread_qc && !params.skip_longread_filtering) {
        if (params.longread_filtering_tool == 'filtlong') {
            tools << "Filtlong (Wick 2019)"
        }
        else if (params.longread_filtering_tool == 'nanoq') {
            tools << "Nanoq (De Coster et al. 2023)"
        }
        else if (params.longread_filtering_tool == 'chopper') {
            tools << "Chopper (De Coster and Rademakers 2023)"
        }
    }
    if (!params.skip_megahit) {
        tools << "MEGAHIT (Li et al. 2015)"
    }
    if (!params.skip_spades) {
        tools << "SPAdes (Bankevich et al. 2012)"
    }
    if (!params.skip_spadeshybrid) {
        tools << "hybridSPAdes (Antipov et al. 2016)"
    }
    if (!params.skip_metamdbg) {
        tools << "metaMDBG (Bourgeade et al. 2024)"
    }
    if (!params.skip_flye) {
        tools << "metaFlye (Kolmogorov et al. 2020)"
    }
    if (!params.skip_quast) {
        tools << "metaQUAST (Mikheenko et al. 2016)"
    }
    if (!params.skip_ale) {
        tools << "ALE (Clark et al. 2013)"
    }
    if (!params.skip_prodigal) {
        tools << "Prodigal (Hyatt et al. 2010)"
    }
    if (!params.skip_prokka) {
        tools << "Prokka (Seemann 2014)"
    }
    if (!params.skip_metaeuk) {
        tools << "MetaEuk (Levy Karin et al. 2020)"
    }
    if (params.run_virus_identification) {
        tools << "geNomad (Camargo et al. 2023)"
    }
    if (!params.skip_binning) {
        if (!params.skip_metabat2) {
            tools << "MetaBAT2 (Kang et al. 2019)"
        }
        if (!params.skip_maxbin2) {
            tools << "MaxBin2 (Wu et al. 2016)"
        }
        if (!params.skip_concoct) {
            tools << "CONCOCT (Alneberg et al. 2014)"
        }
        if (!params.skip_comebin) {
            tools << "COMEBin (Wang et al. 2022)"
        }
        if (!params.skip_metabinner) {
            tools << "MetaBinner (Wang et al. 2023)"
        }
        if (!params.skip_semibin) {
            tools << "SemiBin2 (Pan et al. 2023)"
        }
    }
    if (!params.skip_binqc) {
        if (params.run_busco) {
            tools << "BUSCO (Manni et al. 2021)"
        }
        if (params.run_checkm) {
            tools << "CheckM (Parks et al. 2015)"
        }
        if (params.run_checkm2) {
            tools << "CheckM2 (Chklovski et al. 2023)"
        }
        if (params.refine_bins_dastool) {
            tools << "DAS Tool (Sieber et al. 2018)"
        }
        if (params.run_gunc) {
            tools << "GUNC (Orakov et al. 2021)"
        }
    }
    if (!params.skip_gtdbtk) {
        tools << "GTDB-Tk (Chaumeil et al. 2022)"
    }
    if (params.ancient_dna) {
        tools << "PyDamage (Neukamm et al. 2021)"
        if (!params.skip_ancient_damagecorrection) {
            tools << "FreeBayes (Garrison and Marth 2012)"
            tools << "BCFtools (Danecek et al. 2021)"
        }
    }

    def citation_text = "Tools used in the workflow included: ${tools.join(', ')}."

    return citation_text
}

def toolBibliographyText() {
    def references = [
        "<li>Andrews S. (2010). FastQC. URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</li>",
        "<li>Ewels P, Magnusson M, Lundin S, Kaller M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047-3048. https://doi.org/10.1093/bioinformatics/btw354</li>",
    ]

    if (!params.skip_shortread_qc && !params.skip_clipping) {
        if (params.clip_tool == 'fastp') {
            references << "<li>Chen S, Zhou Y, Chen Y, Gu J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890. https://doi.org/10.1093/bioinformatics/bty560</li>"
        }
        else if (params.clip_tool == 'adapterremoval') {
            references << "<li>Schubert M, Lindgreen S, Orlando L. (2016). AdapterRemoval v2: rapid adapter trimming, identification, and read merging. BMC Research Notes, 9, 88. https://doi.org/10.1186/s13104-016-1900-2</li>"
        }
        else if (params.clip_tool == 'trimmomatic') {
            references << "<li>Bolger AM, Lohse M, Usadel B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics, 30(15), 2114-2120. https://doi.org/10.1093/bioinformatics/btu170</li>"
        }
    }
    if (!params.skip_longread_qc && !params.skip_adapter_trimming) {
        if (params.longread_adaptertrimming_tool == 'porechop') {
            references << "<li>Wick RR. (2017). Porechop. URL: https://github.com/rrwick/Porechop</li>"
        }
        else if (params.longread_adaptertrimming_tool == 'porechop_abi') {
            references << "<li>Wick RR. (2017). Porechop. URL: https://github.com/rrwick/Porechop</li>"
        }
    }
    if (!params.skip_longread_qc && !params.skip_longread_filtering) {
        if (params.longread_filtering_tool == 'filtlong') {
            references << "<li>Wick RR. (2019). Filtlong. URL: https://github.com/rrwick/Filtlong</li>"
        }
        else if (params.longread_filtering_tool == 'nanoq') {
            references << "<li>De Coster W, Rademakers R. (2023). NanoQ: a python package for ultra-fast quality control of Oxford Nanopore sequencing data. Bioinformatics, 39(5), btad311. https://doi.org/10.1093/bioinformatics/btad311</li>"
        }
        else if (params.longread_filtering_tool == 'chopper') {
            references << "<li>De Coster W, Rademakers R. (2023). Chopper: a rust-based quality control tool for nanopore sequencing reads. Journal of Open Source Software, 8(84), 4991. https://doi.org/10.21105/joss.04991</li>"
        }
    }
    if (!params.skip_megahit) {
        references << "<li>Li D, Liu CM, Luo R, Sadakane K, Lam TW. (2015). MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics, 31(10), 1674-1676. https://doi.org/10.1093/bioinformatics/btv033</li>"
    }
    if (!params.skip_spades) {
        references << "<li>Bankevich A, Nurk S, Antipov D, et al. (2012). SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. Journal of Computational Biology, 19(5), 455-477. https://doi.org/10.1089/cmb.2012.0021</li>"
    }
    if (!params.skip_spadeshybrid) {
        references << "<li>Antipov D, Korobeynikov A, McLean JS, Pevzner PA. (2016). hybridSPAdes: an algorithm for hybrid assembly of short and long reads. Bioinformatics, 32(7), 1009-1015. https://doi.org/10.1093/bioinformatics/btv688</li>"
    }
    if (!params.skip_metamdbg) {
        references << "<li>Bourgeade P, Belser C, Bertrand D, et al. (2024). metaMDBG: a scalable long-read metagenome assembler. Nature Biotechnology. https://doi.org/10.1038/s41587-024-02298-5</li>"
    }
    if (!params.skip_flye) {
        references << "<li>Kolmogorov M, Bickhart DM, Behsaz B, et al. (2020). metaFlye: scalable long-read metagenome assembly using repeat graphs. Nature Methods, 17(11), 1103-1110. https://doi.org/10.1038/s41592-020-00971-x</li>"
    }
    if (!params.skip_quast) {
        references << "<li>Mikheenko A, Saveliev V, Gurevich A. (2016). MetaQUAST: evaluation of metagenome assemblies. Bioinformatics, 32(7), 1088-1090. https://doi.org/10.1093/bioinformatics/btv697</li>"
    }
    if (!params.skip_ale) {
        references << "<li>Clark SC, Egan R, Frazier PI, Wang Z. (2013). ALE: a generic assembly likelihood evaluation framework for assessing the accuracy of genome and metagenome assemblies. Bioinformatics, 29(4), 435-443. https://doi.org/10.1093/bioinformatics/bts723</li>"
    }
    if (!params.skip_prodigal) {
        references << "<li>Hyatt D, Chen GL, Locascio PF, Land ML, Larimer FW, Hauser LJ. (2010). Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics, 11, 119. https://doi.org/10.1186/1471-2105-11-119</li>"
    }
    if (!params.skip_prokka) {
        references << "<li>Seemann T. (2014). Prokka: rapid prokaryotic genome annotation. Bioinformatics, 30(14), 2068-2069. https://doi.org/10.1093/bioinformatics/btu153</li>"
    }
    if (!params.skip_metaeuk) {
        references << "<li>Levy Karin E, Mirdita M, Söding J. (2020). MetaEuk-sensitive, high-throughput gene discovery, and annotation for large-scale eukaryotic metagenomics. Microbiome, 8, 48. https://doi.org/10.1186/s40168-020-00808-x</li>"
    }
    if (params.run_virus_identification) {
        references << "<li>Camargo AP, Roux S, Schulz F, et al. (2023). Identification of mobile genetic elements with geNomad. Nature Biotechnology, 41(8), 1080-1089. https://doi.org/10.1038/s41587-023-01953-y</li>"
    }
    if (!params.skip_binning) {
        if (!params.skip_metabat2) {
            references << "<li>Kang DD, Li F, Kirton E, Thomas A, Egan R, An H, Wang Z. (2019). MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. PeerJ, 7, e7359. https://doi.org/10.7717/peerj.7359</li>"
        }
        if (!params.skip_maxbin2) {
            references << "<li>Wu YW, Simmons BA, Singer SW. (2016). MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. Bioinformatics, 32(4), 605-607. https://doi.org/10.1093/bioinformatics/btv638</li>"
        }
        if (!params.skip_concoct) {
            references << "<li>Alneberg J, Bjarnason BS, de Bruijn I, et al. (2014). Binning metagenomic contigs by coverage and composition. Nature Methods, 11(11), 1144-1146. https://doi.org/10.1038/nmeth.3103</li>"
        }
        if (!params.skip_comebin) {
            references << "<li>Wang Z, Niu X, Wang W, et al. (2022). COMEBin: a coverage and composition-based method for metagenomic binning. Research, 2022, 9873837. https://doi.org/10.34133/2022/9873837</li>"
        }
        if (!params.skip_metabinner) {
            references << "<li>Wang Z, Niu X, Zheng Y, et al. (2023). MetaBinner: a high-performance and stand-alone ensemble binning method to recover individual genomes from complex microbial communities. ISME Communications, 3, 7. https://doi.org/10.1038/s43705-022-00235-0</li>"
        }
        if (!params.skip_semibin) {
            references << "<li>Pan S, Zhu C, Zhao XM, Coelho LP. (2023). SemiBin2: self-supervised contrastive learning leads to better MAGs for short- and long-read sequencing. Bioinformatics, 39(Supplement_1), i21-i29. https://doi.org/10.1093/bioinformatics/btad209</li>"
        }
    }
    if (!params.skip_binqc) {
        if (params.run_busco) {
            references << "<li>Manni M, Berkeley MR, Seppey M, Simão FA, Zdobnov EM. (2021). BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes. Molecular Biology and Evolution, 38(10), 4647-4654. https://doi.org/10.1093/molbev/msab199</li>"
        }
        if (params.run_checkm) {
            references << "<li>Parks DH, Imelfort M, Skennerton CT, Hugenholtz P, Tyson GW. (2015). CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Research, 25(7), 1043-1055. https://doi.org/10.1101/gr.186072.114</li>"
        }
        if (params.run_checkm2) {
            references << "<li>Chklovski A, Parks DH, Woodcroft BJ, Tyson GW. (2023). CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. Nature Methods, 20(8), 1203-1212. https://doi.org/10.1038/s41592-023-01940-w</li>"
        }
        if (params.refine_bins_dastool) {
            references << "<li>Sieber CMK, Probst AJ, Sharrar A, et al. (2018). Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. Nature Microbiology, 3(7), 836-843. https://doi.org/10.1038/s41564-018-0171-1</li>"
        }
        if (params.run_gunc) {
            references << "<li>Orakov AN, Fullam A, Coelho LP, et al. (2021). GUNC: detection of chimerism and contamination in prokaryotic genomes. Genome Biology, 22, 178. https://doi.org/10.1186/s13059-021-02393-0</li>"
        }
    }
    if (!params.skip_gtdbtk) {
        references << "<li>Chaumeil PA, Mussig AJ, Hugenholtz P, Parks DH. (2022). GTDB-Tk v2: memory friendly classification with the Genome Taxonomy Database. Bioinformatics, 38(23), 5315-5316. https://doi.org/10.1093/bioinformatics/btac672</li>"
    }
    if (params.ancient_dna) {
        references << "<li>Borry M, Hubner A, Rohrlach AB, Warinner C. (2021). PyDamage: automated ancient damage identification and estimation for contigs in ancient DNA de novo assembly. PeerJ, 9, e11845. https://doi.org/10.7717/peerj.11845</li>"
        if (!params.skip_ancient_damagecorrection) {
            references << "<li>Garrison E, Marth G. (2012). Haplotype-based variant detection from short-read sequencing. arXiv. https://arxiv.org/abs/1207.3907</li>"
            references << "<li>Danecek P, Bonfield JK, Liddle J, et al. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008</li>"
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
