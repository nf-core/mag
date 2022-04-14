//
// This file holds several functions specific to the workflow/mag.nf in the nf-core/mag pipeline
//

class WorkflowMag {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log, hybrid) {
        // Check if binning mapping mode is valid
        if (!['all','group','own'].contains(params.binning_map_mode)) {
            log.error "Invalid parameter '--binning_map_mode ${params.binning_map_mode}'. Valid values are 'all', 'group' or 'own'."
            System.exit(1)
        }
        if (params.coassemble_group && params.binning_map_mode == 'own') {
            log.error "Invalid combination of parameter '--binning_map_mode own' and parameter '--coassemble_group'. Select either 'all' or 'group' mapping mode when performing group-wise co-assembly."
            System.exit(1)
        }

        // Check if specified cpus for SPAdes are available
        if ( params.spades_fix_cpus > params.max_cpus ) {
            log.error "Invalid parameter '--spades_fix_cpus ${params.spades_fix_cpus}', max cpus are '${params.max_cpus}'."
            System.exit(1)
        }
        if ( params.spadeshybrid_fix_cpus > params.max_cpus ) {
            log.error "Invalid parameter '--spadeshybrid_fix_cpus ${params.spadeshybrid_fix_cpus}', max cpus are '${params.max_cpus}'."
            System.exit(1)
        }
        // Check if settings concerning reproducibility of used tools are consistent and print warning if not
        if (params.megahit_fix_cpu_1 || params.spades_fix_cpus != -1 || params.spadeshybrid_fix_cpus != -1){
            if (!params.skip_spades && params.spades_fix_cpus == -1)
                log.warn "At least one assembly process is run with a parameter to ensure reproducible results, but SPAdes not. Consider using the parameter '--spades_fix_cpus'."
            if (hybrid && params.skip_spadeshybrid && params.spadeshybrid_fix_cpus == -1)
                log.warn "At least one assembly process is run with a parameter to ensure reproducible results, but SPAdes hybrid not. Consider using the parameter '--spadeshybrid_fix_cpus'."
            if (!params.skip_megahit && !params.megahit_fix_cpu_1)
                log.warn "At least one assembly process is run with a parameter to ensure reproducible results, but MEGAHIT not. Consider using the parameter '--megahit_fix_cpu_1'."
            if (!params.skip_binning && params.metabat_rng_seed == 0)
                log.warn "At least one assembly process is run with a parameter to ensure reproducible results, but for MetaBAT2 a random seed is specified ('--metabat_rng_seed 0'). Consider specifying a positive seed instead."
        }

        // Check if SPAdes and single_end
        if ( (!params.skip_spades || !params.skip_spadeshybrid) && params.single_end) {
            log.warn "metaSPAdes does not support single-end data. SPAdes will be skipped."
        }

        // Check if parameters for host contamination removal are valid
        if ( params.host_fasta && params.host_genome) {
            log.error "Both host fasta reference and iGenomes genome are specified to remove host contamination! Invalid combination, please specify either --host_fasta or --host_genome."
            System.exit(1)
        }
        if ( hybrid && (params.host_fasta || params.host_genome) ) {
            log.warn "Host read removal is only applied to short reads. Long reads might be filtered indirectly by Filtlong, which is set to use read qualities estimated based on k-mer matches to the short, already filtered reads."
            if ( params.longreads_length_weight > 1 ) {
                log.warn "The parameter --longreads_length_weight is ${params.longreads_length_weight}, causing the read length being more important for long read filtering than the read quality. Set --longreads_length_weight to 1 in order to assign equal weights."
            }
        }
        if ( params.host_genome ) {
            if (!params.genomes) {
                log.error "No config file containing genomes provided!"
                System.exit(1)
            }
            // Check if host genome exists in the config file
            if (!params.genomes.containsKey(params.host_genome)) {
                log.error "=============================================================================\n" +
                        "  Host genome '${params.host_genome}' not found in any config files provided to the pipeline.\n" +
                        "  Currently, the available genome keys are:\n" +
                        "  ${params.genomes.keySet().join(", ")}\n" +
                        "==================================================================================="
                System.exit(1)
            }
            if ( !params.genomes[params.host_genome].fasta ) {
                log.error "No fasta file specified for the host genome ${params.host_genome}!"
                System.exit(1)
            }
            if ( !params.genomes[params.host_genome].bowtie2 ) {
                log.error "No Bowtie 2 index file specified for the host genome ${params.host_genome}!"
                System.exit(1)
            }
        }

        // Check if t least two binners where applied in order to run DAS Tool for bin refinment
        // (needs to be adjusted in case additional binners are added)
        if (params.refine_bins_dastool && params.skip_metabat2 ){
            log.error "Both --refine_bins_dastool and --skip_metabat2 are specified! Invalid combination, bin refinement requires MetaBAT2 and MaxBin2 binning results."
            System.exit(1)
        }
        if (params.refine_bins_dastool && params.skip_maxbin2 ){
            log.error "Both --refine_bins_dastool and --skip_maxbin2 are specified! Invalid combination, bin refinement requires MetaBAT2 and MaxBin2 binning results."
            System.exit(1)
        }

        // Check if BUSCO parameters combinations are valid
        if (params.skip_busco){
            if (params.busco_reference) {
                log.error "Both --skip_busco and --busco_reference are specified! Invalid combination, please specify either --skip_busco or --busco_reference."
                System.exit(1)
            }
            if (params.busco_download_path) {
                log.error "Both --skip_busco and --busco_download_path are specified! Invalid combination, please specify either --skip_busco or --busco_download_path."
                System.exit(1)
            }
            if (params.busco_auto_lineage_prok) {
                log.error "Both --skip_busco and --busco_auto_lineage_prok are specified! Invalid combination, please specify either --skip_busco or --busco_auto_lineage_prok."
                System.exit(1)
            }
        }
        if (params.busco_reference && params.busco_download_path) {
            log.error "Both --busco_reference and --busco_download_path are specified! Invalid combination, please specify either --busco_reference or --busco_download_path."
            System.exit(1)
        }
        if (params.busco_auto_lineage_prok && params.busco_reference) {
            log.error "Both --busco_auto_lineage_prok and --busco_reference are specified! Invalid combination, please specify either --busco_auto_lineage_prok or --busco_reference."
            System.exit(1)
        }

        if (params.skip_busco && params.gtdb) {
            log.warn "--skip_busco and --gtdb are specified! GTDB-tk will be omitted because GTDB-tk bin classification requires bin filtering based on BUSCO QC results to avoid GTDB-tk errors."
        }

        // Check if CAT parameters are valid
        if (params.cat_db && params.cat_db_generate) {
            log.error "Invalid combination of parameters --cat_db and --cat_db_generate is specified! Please specify either --cat_db or --cat_db_generate."
            System.exit(1)
        }
        if (params.save_cat_db && !params.cat_db_generate) {
            log.error "Invalid parameter combination: parameter --save_cat_db specified, but not --cat_db_generate! Note also that the parameter --save_cat_db does not work in combination with --cat_db."
            System.exit(1)
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
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }
}
