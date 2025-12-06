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
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    _monochrome_logs  // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

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
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

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
    if (params.refine_bins_dastool && !([params.skip_metabat2, params.skip_maxbin2, params.skip_concoct].count(false) > 1)) {
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
            log.warn('[nf-core/mag]: --skip_binqc is specified, but --skip_gtdbtk is explictly set to run! GTDB-tk will be omitted because GTDB-tk bin classification requires bin filtering based on BUSCO or CheckM QC results to avoid GTDB-tk errors.')
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

    // Check MetaEuk db paramaters
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
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(meta, sr1, sr2, lr) {

    if ((!sr2 && !lr) && !params.single_end) {
        error("[nf-core/mag] ERROR: Single-end data must be executed with `--single_end`. Note that it is not possible to mix single- and paired-end data in one run! Check input TSV for sample: ${meta.id}")
    }
    if (sr2 && params.single_end) {
        error("[nf-core/mag] ERROR: Paired-end data must be executed without `--single_end`. Note that it is not possible to mix single- and paired-end data in one run! Check input TSV for sample: ${meta.id}")
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
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
        "Tools used in the workflow included:",
        "FastQC (Andrews 2010),",
        "MultiQC (Ewels et al. 2016)",
        ".",
    ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
        "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
        "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>",
    ].join(' ').trim()

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
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine = new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
