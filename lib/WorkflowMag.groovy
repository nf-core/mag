//
// This file holds several functions specific to the workflow/mag.nf in the nf-core/mag pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowMag {

    //
    // Check and validate parameters
    //

    public static void initialise(params, log, hybrid) {
        // Check if binning mapping mode is valid
        if (!['all', 'group', 'own'].contains(params.binning_map_mode)) {
            Nextflow.error("Invalid parameter '--binning_map_mode ${params.binning_map_mode}'. Valid values are 'all', 'group' or 'own'.")
        }
        if (params.coassemble_group && params.binning_map_mode == 'own') {
            Nextflow.error("Invalid combination of parameter '--binning_map_mode own' and parameter '--coassemble_group'. Select either 'all' or 'group' mapping mode when performing group-wise co-assembly.")
        }
        if (params.ancient_dna && params.binning_map_mode != 'own') {
            Nextflow.error("Invalid combination of parameter '--binning_map_mode' and parameter '--ancient_dna'. Ancient DNA mode can only be executed with --binning_map_mode own. You supplied: --binning_map_mode ${params.binning_map_mode}")
        }

        // Check if specified cpus for SPAdes are available
        if ( params.spades_fix_cpus > params.max_cpus ) {
            Nextflow.error("Invalid parameter '--spades_fix_cpus ${params.spades_fix_cpus}', max cpus are '${params.max_cpus}'.")
        }
        if ( params.spadeshybrid_fix_cpus > params.max_cpus ) {
            Nextflow.error("Invalid parameter '--spadeshybrid_fix_cpus ${params.spadeshybrid_fix_cpus}', max cpus are '${params.max_cpus}'.")
        }
        // Check if settings concerning reproducibility of used tools are consistent and print warning if not
        if (params.megahit_fix_cpu_1 || params.spades_fix_cpus != -1 || params.spadeshybrid_fix_cpus != -1) {
            if (!params.skip_spades && params.spades_fix_cpus == -1) {
                log.warn "At least one assembly process is run with a parameter to ensure reproducible results, but SPAdes not. Consider using the parameter '--spades_fix_cpus'."
            }
            if (hybrid && params.skip_spadeshybrid && params.spadeshybrid_fix_cpus == -1) {
                log.warn "At least one assembly process is run with a parameter to ensure reproducible results, but SPAdes hybrid not. Consider using the parameter '--spadeshybrid_fix_cpus'."
            }
            if (!params.skip_megahit && !params.megahit_fix_cpu_1) {
                log.warn "At least one assembly process is run with a parameter to ensure reproducible results, but MEGAHIT not. Consider using the parameter '--megahit_fix_cpu_1'."
            }
            if (!params.skip_binning && params.metabat_rng_seed == 0) {
                log.warn "At least one assembly process is run with a parameter to ensure reproducible results, but for MetaBAT2 a random seed is specified ('--metabat_rng_seed 0'). Consider specifying a positive seed instead."
            }
        }

        // Check if SPAdes and single_end
        if ( (!params.skip_spades || !params.skip_spadeshybrid) && params.single_end) {
            log.warn 'metaSPAdes does not support single-end data. SPAdes will be skipped.'
        }

        // Check if parameters for host contamination removal are valid
        if ( params.host_fasta && params.host_genome) {
            Nextflow.error('Both host fasta reference and iGenomes genome are specified to remove host contamination! Invalid combination, please specify either --host_fasta or --host_genome.')
        }
        if ( hybrid && (params.host_fasta || params.host_genome) ) {
            log.warn 'Host read removal is only applied to short reads. Long reads might be filtered indirectly by Filtlong, which is set to use read qualities estimated based on k-mer matches to the short, already filtered reads.'
            if ( params.longreads_length_weight > 1 ) {
                log.warn "The parameter --longreads_length_weight is ${params.longreads_length_weight}, causing the read length being more important for long read filtering than the read quality. Set --longreads_length_weight to 1 in order to assign equal weights."
            }
        }
        if ( params.host_genome ) {
            if (!params.genomes) {
                Nextflow.error('No config file containing genomes provided!')
            }
            // Check if host genome exists in the config file
            if (!params.genomes.containsKey(params.host_genome)) {
                Nextflow.error('=============================================================================\n' +
                        "  Host genome '${params.host_genome}' not found in any config files provided to the pipeline.\n" +
                        '  Currently, the available genome keys are:\n' +
                        "  ${params.genomes.keySet().join(', ')}\n" +
                        '===================================================================================')
            }
            if ( !params.genomes[params.host_genome].fasta ) {
                Nextflow.error("No fasta file specified for the host genome ${params.host_genome}!")
            }
            if ( !params.genomes[params.host_genome].bowtie2 ) {
                Nextflow.error("No Bowtie 2 index file specified for the host genome ${params.host_genome}!")
            }
        }

        // Check MetaBAT2 inputs
        if ( !params.skip_metabat2 && params.min_contig_size < 1500 ) {
            log.warn "Specified min. contig size under minimum for MetaBAT2. MetaBAT2 will be run with 1500 (other binners not affected). You supplied: --min_contig_size ${params.min_contig_size}"
        }

        // Check more than one binner is run for bin refinement  (required DAS by Tool)
        // If the number of run binners (i.e., number of not-skipped) is more than one, otherwise throw an error
        if ( params.refine_bins_dastool && !([ params.skip_metabat2, params.skip_maxbin2, params.skip_concoct ].count(false) > 1) ) {
            Nextflow.error('Bin refinement with --refine_bins_dastool requires at least two binners to be running (not skipped). Check input.')
        }

        // Check that bin refinement is actually turned on if any of the refined bins are requested for downstream
        if (!params.refine_bins_dastool && params.postbinning_input != 'raw_bins_only') {
            Nextflow.error("The parameter '--postbinning_input ${ params.postbinning_input }' for downstream steps can only be specified if bin refinement is activated with --refine_bins_dastool! Check input.")
        }

        // Check if BUSCO parameters combinations are valid
        if (params.skip_binqc && params.binqc_tool == 'checkm') {
            Nextflow.error('Both --skip_binqc and --binqc_tool \'checkm\' are specified! Invalid combination, please specify either --skip_binqc or --binqc_tool.')
        }
        if (params.skip_binqc) {
            if (params.busco_db) {
                Nextflow.error('Both --skip_binqc and --busco_db are specified! Invalid combination, please specify either --skip_binqc or --binqc_tool \'busco\' with --busco_db.')
            }
            if (params.busco_auto_lineage_prok) {
                Nextflow.error('Both --skip_binqc and --busco_auto_lineage_prok are specified! Invalid combination, please specify either --skip_binqc or --binqc_tool \'busco\' with --busco_auto_lineage_prok.')
            }
        }

        if (params.skip_binqc && !params.skip_gtdbtk) {
            log.warn '--skip_binqc is specified, but --skip_gtdbtk is explictly set to run! GTDB-tk will be omitted because GTDB-tk bin classification requires bin filtering based on BUSCO or CheckM QC results to avoid GTDB-tk errors.'
        }

        // Check if CAT parameters are valid
        if (params.cat_db && params.cat_db_generate) {
            Nextflow.error('Invalid combination of parameters --cat_db and --cat_db_generate is specified! Please specify either --cat_db or --cat_db_generate.')
        }
        if (params.save_cat_db && !params.cat_db_generate) {
            Nextflow.error('Invalid parameter combination: parameter --save_cat_db specified, but not --cat_db_generate! Note also that the parameter --save_cat_db does not work in combination with --cat_db.')
        }

        // Chech MetaEuk db paramaters
        if (params.metaeuk_mmseqs_db && params.metaeuk_db) {
            Nextflow.error('Invalid parameter combination: both --metaeuk_mmseqs_db and --metaeuk_db are specified! Please specify either --metaeuk_mmseqs_db or --metaeuk_db.')
        }
        if (params.save_mmseqs_db && !params.metaeuk_mmseqs_db) {
            Nextflow.error('Invalid parameter combination: --save_mmseqs_db supplied but no database has been requested for download with --metaeuk_mmseqs_db!')
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += '    </dl>\n'
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/', '-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += 'data: |\n'
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Generate methods description for MultiQC
    //

    public static String toolCitationText(params) {

        // TODO Optionally add in-text citation tools to this list.
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

    public static String toolBibliographyText(params) {

        // TODO Optionally add bibliographic entries to this list.
        // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
        // Uncomment function in methodsDescriptionText to render in MultiQC report
        def reference_text = [
                "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
                "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
            ].join(' ').trim()

        return reference_text
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml, params) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta['manifest_map'] = run_workflow.manifest.toMap()

        // Pipeline DOI
        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        // Tool references
        meta["tool_citations"] = ""
        meta["tool_bibliography"] = ""

        // TODO Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
        //meta["tool_citations"] = toolCitationText(params).replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
        //meta["tool_bibliography"] = toolBibliographyText(params)


        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.error(error_string)
        }
    }

}
