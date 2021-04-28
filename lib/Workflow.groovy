/*
 * This file holds several functions specific to the pipeline.
 */

class Workflow {

    // Citation string
    private static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
               "* The pipeline\n" + 
               "  https://doi.org/10.5281/zenodo.3589527\n\n" +
               "* The nf-core framework\n" +
               "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
               "* Software dependencies\n" +
               "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }

    static void validateMainParams(workflow, params, json_schema, log) {
        if (params.validate_params) {
            NfcoreSchema.validateParameters(params, json_schema, log)
        }

        // Check that conda channels are set-up correctly
        if (params.enable_conda) {
            Checks.checkCondaChannels(log)
        }

        // Check AWS batch settings
        Checks.awsBatch(workflow, params)

        // Check the hostnames against configured profiles
        Checks.hostName(workflow, params, log)
    }

    static void validateWorkflowParams(params, log, hybrid) {
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
            log.error "Both host fasta reference and iGenomes genome are specififed to remove host contamination! Invalid combination, please specify either --host_fasta or --host_genome."
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

        // Check if BUSCO parameters combinations are valid
        if (params.skip_busco){
            if (params.busco_reference) {
                log.error "Both --skip_busco and --busco_reference are specififed! Invalid combination, please specify either --skip_busco or --busco_reference."
                System.exit(1)
            }
            if (params.busco_download_path) {
                log.error "Both --skip_busco and --busco_download_path are specififed! Invalid combination, please specify either --skip_busco or --busco_download_path."
                System.exit(1)
            }
            if (params.busco_auto_lineage_prok) {
                log.error "Both --skip_busco and --busco_auto_lineage_prok are specififed! Invalid combination, please specify either --skip_busco or --busco_auto_lineage_prok."
                System.exit(1)
            }
        }
        if (params.busco_reference && params.busco_download_path) {
            log.error "Both --busco_reference and --busco_download_path are specififed! Invalid combination, please specify either --busco_reference or --busco_download_path."
            System.exit(1)
        }
        if (params.busco_auto_lineage_prok && params.busco_reference) {
            log.error "Both --busco_auto_lineage_prok and --busco_reference are specififed! Invalid combination, please specify either --busco_auto_lineage_prok or --busco_reference."
            System.exit(1)
        }
    }

    /*
     * Get workflow summary for MultiQC
     */
    static String paramsSummaryMultiqc(workflow, summary) {
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
