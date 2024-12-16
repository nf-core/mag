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

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE(
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN(
        workflow,
        validate_params,
        null
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
    ch_samplesheet = Channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map {
            validateInputSamplesheet(it[0], it[1], it[2], it[3])
        }

    // Prepare FASTQs channel and separate short and long reads and prepare
    ch_raw_short_reads = ch_samplesheet.map { meta, sr1, sr2, _lr ->
        meta.run = meta.run == [] ? "0" : meta.run
        meta.single_end = params.single_end

        if (params.single_end) {
            return [meta, [sr1]]
        }
        else {
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
    ch_raw_long_reads.map {
        if (it) {
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
        ch_input_assemblies = Channel.fromList(samplesheetToList(params.assembly_input, "${projectDir}/assets/schema_assembly_input.json"))
    }

    // Prepare ASSEMBLY input channel
    if (params.assembly_input) {
        ch_input_assemblies.map { meta, fasta ->
            return [meta + [id: params.coassemble_group ? "group-${meta.group}" : meta.id], [fasta]]
        }
    }
    else {
        ch_input_assemblies = Channel.empty()
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

    // Check if settings concerning reproducibility of used tools are consistent and print warning if not
    if (params.megahit_fix_cpu_1 || params.spades_fix_cpus != -1 || params.spadeshybrid_fix_cpus != -1) {
        if (!params.skip_spades && params.spades_fix_cpus == -1) {
            log.warn("[nf-core/mag]: At least one assembly process is run with a parameter to ensure reproducible results, but SPAdes not. Consider using the parameter '--spades_fix_cpus'.")
        }
        if (hybrid && params.skip_spadeshybrid && params.spadeshybrid_fix_cpus == -1) {
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
    if (hybrid && (params.host_fasta || params.host_genome)) {
        log.warn('[nf-core/mag]: Host read removal is only applied to short reads. Long reads might be filtered indirectly by Filtlong, which is set to use read qualities estimated based on k-mer matches to the short, already filtered reads.')
        if (params.longreads_length_weight > 1) {
            log.warn("[nf-core/mag]: The parameter --longreads_length_weight is ${params.longreads_length_weight}, causing the read length being more important for long read filtering than the read quality. Set --longreads_length_weight to 1 in order to assign equal weights.")
        }
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

    // Check if BUSCO parameters combinations are valid
    if (params.skip_binqc && params.binqc_tool == 'checkm') {
        error("[nf-core/mag] ERROR: Both --skip_binqc and --binqc_tool 'checkm' are specified! Invalid combination, please specify either --skip_binqc or --binqc_tool.")
    }
    if (params.skip_binqc) {
        if (params.busco_db) {
            error("[nf-core/mag] ERROR: Both --skip_binqc and --busco_db are specified! Invalid combination, please specify either --skip_binqc or --binqc_tool 'busco' with --busco_db.")
        }
        if (params.busco_auto_lineage_prok) {
            error("[nf-core/mag] ERROR: Both --skip_binqc and --busco_auto_lineage_prok are specified! Invalid combination, please specify either --skip_binqc or --binqc_tool 'busco' with --busco_auto_lineage_prok.")
        }
    }

    if (params.skip_binqc && !params.skip_gtdbtk) {
        log.warn('[nf-core/mag]: --skip_binqc is specified, but --skip_gtdbtk is explictly set to run! GTDB-tk will be omitted because GTDB-tk bin classification requires bin filtering based on BUSCO or CheckM QC results to avoid GTDB-tk errors.')
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
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(meta, sr1, sr2, lr) {

    if (!sr2 && !params.single_end) {
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
def genomeExistsError() {
    if (params.genomes && params.host_genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" + "  Genome '${params.host_genome}' not found in any config files provided to the pipeline.\n" + "  Currently, the available genome keys are:\n" + "  ${params.host_genomes.keySet().join(", ")}\n" + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
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
        "."
    ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
        "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
        "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
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
