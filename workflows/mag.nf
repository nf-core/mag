/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                                               } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                                      } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                                  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                                } from '../subworkflows/local/utils_nfcore_mag_pipeline'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { BINNING_PREPARATION                                   } from '../subworkflows/local/binning_preparation'
include { BINNING                                               } from '../subworkflows/local/binning'
include { BIN_QC                                                } from '../subworkflows/local/bin_qc'
include { BINNING_REFINEMENT                                    } from '../subworkflows/local/binning_refinement'
include { VIRUS_IDENTIFICATION                                  } from '../subworkflows/local/virus_identification'
include { GTDBTK                                                } from '../subworkflows/local/gtdbtk'
include { ANCIENT_DNA_ASSEMBLY_VALIDATION                       } from '../subworkflows/local/ancient_dna'
include { DOMAIN_CLASSIFICATION                                 } from '../subworkflows/local/domain_classification'
include { DEPTHS                                                } from '../subworkflows/local/depths'
include { LONGREAD_PREPROCESSING                                } from '../subworkflows/local/longread_preprocessing'
include { SHORTREAD_PREPROCESSING                               } from '../subworkflows/local/shortread_preprocessing'

//
// MODULE: Installed directly from nf-core/modules
//
include { UNTAR as CENTRIFUGEDB_UNTAR                           } from '../modules/nf-core/untar/main'
include { CENTRIFUGE_CENTRIFUGE                                 } from '../modules/nf-core/centrifuge/centrifuge/main'
include { CENTRIFUGE_KREPORT                                    } from '../modules/nf-core/centrifuge/kreport/main'
include { KRONA_KRONADB                                         } from '../modules/nf-core/krona/kronadb/main'
include { KRONA_KTIMPORTTAXONOMY                                } from '../modules/nf-core/krona/ktimporttaxonomy/main'
include { KRAKENTOOLS_KREPORT2KRONA as KREPORT2KRONA_CENTRIFUGE } from '../modules/nf-core/krakentools/kreport2krona/main'
include { MEGAHIT                                               } from '../modules/nf-core/megahit/main'
include { SPADES as METASPADES                                  } from '../modules/nf-core/spades/main'
include { SPADES as METASPADESHYBRID                            } from '../modules/nf-core/spades/main'
include { GUNZIP as GUNZIP_ASSEMBLIES                           } from '../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_ASSEMBLYINPUT                        } from '../modules/nf-core/gunzip'
include { PRODIGAL                                              } from '../modules/nf-core/prodigal/main'
include { PROKKA                                                } from '../modules/nf-core/prokka/main'
include { MMSEQS_DATABASES                                      } from '../modules/nf-core/mmseqs/databases/main'
include { METAEUK_EASYPREDICT                                   } from '../modules/nf-core/metaeuk/easypredict/main'

//
// MODULE: Local to the pipeline
//
include { KRAKEN2_DB_PREPARATION                                } from '../modules/local/kraken2_db_preparation'
include { KRAKEN2                                               } from '../modules/local/kraken2'
include { POOL_SINGLE_READS as POOL_SHORT_SINGLE_READS          } from '../modules/local/pool_single_reads'
include { POOL_PAIRED_READS                                     } from '../modules/local/pool_paired_reads'
include { POOL_SINGLE_READS as POOL_LONG_READS                  } from '../modules/local/pool_single_reads'
include { QUAST                                                 } from '../modules/local/quast'
include { QUAST_BINS                                            } from '../modules/local/quast_bins'
include { QUAST_BINS_SUMMARY                                    } from '../modules/local/quast_bins_summary'
include { CAT_DB                                                } from '../modules/local/cat_db'
include { CAT_DB_GENERATE                                       } from '../modules/local/cat_db_generate'
include { CAT                                                   } from '../modules/local/cat'
include { CAT_SUMMARY                                           } from '../modules/local/cat_summary'
include { BIN_SUMMARY                                           } from '../modules/local/bin_summary'
include { COMBINE_TSV as COMBINE_SUMMARY_TSV                    } from '../modules/local/combine_tsv'

workflow MAG {
    take:
    ch_raw_short_reads  // channel: samplesheet read in from --input
    ch_raw_long_reads
    ch_input_assemblies

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ////////////////////////////////////////////////////
    /* --  Create channel for reference databases  -- */
    ////////////////////////////////////////////////////

    if (params.host_genome) {
        host_fasta = params.genomes[params.host_genome].fasta ?: false
        ch_host_fasta = Channel.value(file("${host_fasta}"))
        host_bowtie2index = params.genomes[params.host_genome].bowtie2 ?: false
        ch_host_bowtie2index = Channel.fromPath("${host_bowtie2index}", checkIfExists: true).first()
    }
    else if (params.host_fasta) {
        ch_host_fasta = Channel.fromPath("${params.host_fasta}", checkIfExists: true).first() ?: false

        if (params.host_fasta_bowtie2index) {
            ch_host_bowtie2index = Channel.fromPath("${params.host_fasta_bowtie2index}", checkIfExists: true).first()
        }
        else {
            ch_host_bowtie2index = Channel.empty()
        }
    }
    else {
        ch_host_fasta = Channel.empty()
        ch_host_bowtie2index = Channel.empty()
    }

    if (params.kraken2_db) {
        ch_kraken2_db_file = file(params.kraken2_db, checkIfExists: true)
    }
    else {
        ch_kraken2_db_file = []
    }

    if (params.cat_db) {
        ch_cat_db_file = Channel.value(file("${params.cat_db}"))
    }
    else {
        ch_cat_db_file = Channel.empty()
    }

    if (params.krona_db) {
        ch_krona_db_file = Channel.value(file("${params.krona_db}"))
    }
    else {
        ch_krona_db_file = Channel.empty()
    }

    if (!params.keep_phix) {
        ch_phix_db_file = Channel.value(file("${params.phix_reference}"))
    }
    else {
        ch_phix_db_file = Channel.empty()
    }

    if (!params.keep_lambda) {
        ch_lambda_db = Channel.value(file("${params.lambda_reference}"))
    }
    else {
        ch_lambda_db = Channel.empty()
    }

    if (params.genomad_db) {
        ch_genomad_db = file(params.genomad_db, checkIfExists: true)
    }
    else {
        ch_genomad_db = Channel.empty()
    }

    gtdb = params.skip_binqc || params.skip_gtdbtk ? false : params.gtdb_db

    if (gtdb) {
        gtdb = file("${gtdb}", checkIfExists: true)
        gtdb_mash = params.gtdb_mash ? file("${params.gtdb_mash}", checkIfExists: true) : []
    }
    else {
        gtdb = []
    }

    if (params.metaeuk_db && !params.skip_metaeuk) {
        ch_metaeuk_db = Channel.value(file("${params.metaeuk_db}", checkIfExists: true))
    }
    else {
        ch_metaeuk_db = Channel.empty()
    }

    // Get mmseqs db for MetaEuk if requested
    if (!params.skip_metaeuk && params.metaeuk_mmseqs_db) {
        MMSEQS_DATABASES(params.metaeuk_mmseqs_db)
        ch_metaeuk_db = MMSEQS_DATABASES.out.database
        ch_versions = ch_versions.mix(MMSEQS_DATABASES.out.versions)
    }

    /*
    ================================================================================
                                    Preprocessing and QC for short reads
    ================================================================================
    */

    if (!params.assembly_input) {
        SHORTREAD_PREPROCESSING(
            ch_raw_short_reads,
            ch_host_fasta,
            ch_host_bowtie2index,
            ch_phix_db_file,
        )

        ch_versions = ch_versions.mix(SHORTREAD_PREPROCESSING.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(SHORTREAD_PREPROCESSING.out.multiqc_files.collect { it[1] }.ifEmpty([]))
        ch_short_reads = SHORTREAD_PREPROCESSING.out.short_reads
        ch_short_reads_assembly = SHORTREAD_PREPROCESSING.out.short_reads_assembly
    }
    else {
        ch_short_reads = ch_raw_short_reads.map { meta, reads ->
            def meta_new = meta - meta.subMap('run')
            [meta_new, reads]
        }
    }

    /*
    ================================================================================
                                    Preprocessing and QC for long reads
    ================================================================================
    */

    LONGREAD_PREPROCESSING(
        ch_raw_long_reads,
        ch_short_reads,
        ch_lambda_db,
    )

    ch_versions = ch_versions.mix(LONGREAD_PREPROCESSING.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(LONGREAD_PREPROCESSING.out.multiqc_files.collect { it[1] }.ifEmpty([]))
    ch_long_reads = LONGREAD_PREPROCESSING.out.long_reads

    /*
    ================================================================================
                                    Taxonomic information
    ================================================================================
    */

    // Centrifuge
    if (!params.centrifuge_db) {
        ch_db_for_centrifuge = Channel.empty()
    }
    else {
        if (file(params.centrifuge_db).isDirectory()) {
            ch_db_for_centrifuge = Channel.of(file(params.centrifuge_db, checkIfExists: true))
        }
        else {
            ch_db_for_centrifuge = CENTRIFUGEDB_UNTAR(Channel.of([[id: 'db'], file(params.centrifuge_db, checkIfExists: true)])).untar.map { it[1] }.first()
            ch_versions = ch_versions.mix(CENTRIFUGEDB_UNTAR.out.versions.first())
        }
    }

    CENTRIFUGE_CENTRIFUGE(
        ch_short_reads,
        ch_db_for_centrifuge,
        false,
        false,
    )
    ch_versions = ch_versions.mix(CENTRIFUGE_CENTRIFUGE.out.versions.first())

    CENTRIFUGE_KREPORT(CENTRIFUGE_CENTRIFUGE.out.results, ch_db_for_centrifuge)
    ch_versions = ch_versions.mix(CENTRIFUGE_KREPORT.out.versions.first())

    // Kraken2
    if (!ch_kraken2_db_file.isEmpty()) {
        if (ch_kraken2_db_file.extension in ['gz', 'tgz']) {
            // Expects to be tar.gz!
            ch_db_for_kraken2 = KRAKEN2_DB_PREPARATION(ch_kraken2_db_file).db
        }
        else if (ch_kraken2_db_file.isDirectory()) {
            ch_db_for_kraken2 = Channel
                .fromPath("${ch_kraken2_db_file}/*.k2d")
                .collect()
                .map { file ->
                    if (file.size() >= 3) {
                        def db_name = file[0].getParent().getName()
                        [db_name, file]
                    }
                    else {
                        error("Kraken2 requires '{hash,opts,taxo}.k2d' files.")
                    }
                }
        }
        else {
            ch_db_for_kraken2 = Channel.empty()
        }
    }
    else {
        ch_db_for_kraken2 = Channel.empty()
    }

    KRAKEN2(
        ch_short_reads,
        ch_db_for_kraken2,
    )
    ch_versions = ch_versions.mix(KRAKEN2.out.versions.first())

    if ((params.centrifuge_db || params.kraken2_db) && !params.skip_krona) {
        if (params.krona_db) {
            ch_krona_db = ch_krona_db_file
        }
        else {
            KRONA_KRONADB()
            ch_krona_db = KRONA_KRONADB.out.db
            ch_versions = ch_versions.mix(KRONA_KRONADB.out.versions)
        }

        if (params.centrifuge_db) {
            ch_centrifuge_for_krona = KREPORT2KRONA_CENTRIFUGE(CENTRIFUGE_KREPORT.out.kreport).txt.map { meta, files -> ['centrifuge', meta, files] }
            ch_versions = ch_versions.mix(KREPORT2KRONA_CENTRIFUGE.out.versions.first())
        }
        else {
            ch_centrifuge_for_krona = Channel.empty()
        }

        // Join together for Krona
        ch_tax_classifications = ch_centrifuge_for_krona
            .mix(KRAKEN2.out.results_for_krona)
            .map { classifier, meta, report ->
                def meta_new = meta + [classifier: classifier]
                [meta_new, report]
            }

        KRONA_KTIMPORTTAXONOMY(
            ch_tax_classifications,
            ch_krona_db,
        )
        ch_versions = ch_versions.mix(KRONA_KTIMPORTTAXONOMY.out.versions.first())
    }

    /*
    ================================================================================
                                    Assembly
    ================================================================================
    */

    if (!params.assembly_input) {

        // Co-assembly preparation: grouping for MEGAHIT and for pooling for SPAdes
        if (params.coassemble_group) {
            // short reads
            // group and set group as new id
            ch_short_reads_grouped = ch_short_reads_assembly
                .map { meta, reads -> [meta.group, meta, reads] }
                .groupTuple(by: 0)
                .map { group, _metas, reads ->
                    def assemble_as_single = params.single_end || (params.bbnorm && params.coassemble_group)
                    def meta = [:]
                    meta.id = "group-${group}"
                    meta.group = group
                    meta.single_end = assemble_as_single
                    if (assemble_as_single) {
                        [meta, reads.collect { it }, []]
                    }
                    else {
                        [meta, reads.collect { it[0] }, reads.collect { it[1] }]
                    }
                }
            // long reads
            // group and set group as new id
            ch_long_reads_grouped = ch_long_reads
                .map { meta, reads -> [meta.group, meta, reads] }
                .groupTuple(by: 0)
                .map { group, _metas, reads ->
                    def meta = [:]
                    meta.id = "group-${group}"
                    meta.group = group
                    [meta, reads.collect { it }]
                }
        }
        else {
            ch_short_reads_grouped = ch_short_reads_assembly
                .filter { it[0].single_end }
                .map { meta, reads -> [meta, [reads], []] }
                .mix(
                    ch_short_reads_assembly.filter { !it[0].single_end }.map { meta, reads -> [meta, [reads[0]], [reads[1]]] }
                )
            ch_long_reads_grouped = ch_long_reads
        }

        if (!params.skip_spades || !params.skip_spadeshybrid) {
            if (params.coassemble_group) {
                if (params.bbnorm) {
                    ch_short_reads_spades = ch_short_reads_grouped.map { [it[0], it[1]] }
                }
                else {
                    POOL_SHORT_SINGLE_READS(
                        ch_short_reads_grouped.filter { it[0].single_end }
                    )
                    POOL_PAIRED_READS(
                        ch_short_reads_grouped.filter { !it[0].single_end }
                    )
                    ch_short_reads_spades = POOL_SHORT_SINGLE_READS.out.reads.mix(POOL_PAIRED_READS.out.reads)
                }
            }
            else {
                ch_short_reads_spades = ch_short_reads_assembly
            }
            // long reads
            if (!params.single_end && !params.skip_spadeshybrid) {
                POOL_LONG_READS(ch_long_reads_grouped)
                ch_long_reads_spades = POOL_LONG_READS.out.reads
            }
            else {
                ch_long_reads_spades = Channel.empty()
            }
        }
        else {
            ch_short_reads_spades = Channel.empty()
            ch_long_reads_spades = Channel.empty()
        }

        // Assembly

        ch_assembled_contigs = Channel.empty()

        if (!params.single_end && !params.skip_spades) {
            METASPADES(ch_short_reads_spades.map { meta, reads -> [meta, reads, [], []] }, [], [])
            ch_spades_assemblies = (params.spades_downstreaminput == 'contigs' ? METASPADES.out.contigs : METASPADES.out.scaffolds).map { meta, assembly ->
                def meta_new = meta + [assembler: 'SPAdes']
                [meta_new, assembly]
            }
            ch_assembled_contigs = ch_assembled_contigs.mix(ch_spades_assemblies)
            ch_versions = ch_versions.mix(METASPADES.out.versions.first())
        }

        if (!params.single_end && !params.skip_spadeshybrid) {
            ch_short_reads_spades_tmp = ch_short_reads_spades.map { meta, reads -> [meta.id, meta, reads] }

            ch_reads_spadeshybrid = ch_long_reads_spades
                .map { meta, reads -> [meta.id, meta, reads] }
                .combine(ch_short_reads_spades_tmp, by: 0)
                .map { _id, _meta_long, long_reads, meta_short, short_reads -> [meta_short, short_reads, [], long_reads] }

            METASPADESHYBRID(ch_reads_spadeshybrid, [], [])
            ch_spadeshybrid_assemblies = METASPADESHYBRID.out.scaffolds.map { meta, assembly ->
                def meta_new = meta + [assembler: "SPAdesHybrid"]
                [meta_new, assembly]
            }
            ch_assembled_contigs = ch_assembled_contigs.mix(ch_spadeshybrid_assemblies)
            ch_versions = ch_versions.mix(METASPADESHYBRID.out.versions.first())
        }

        if (!params.skip_megahit) {
            MEGAHIT(ch_short_reads_grouped)
            ch_megahit_assemblies = MEGAHIT.out.contigs.map { meta, assembly ->
                def meta_new = meta + [assembler: 'MEGAHIT']
                [meta_new, assembly]
            }
            ch_assembled_contigs = ch_assembled_contigs.mix(ch_megahit_assemblies)
            ch_versions = ch_versions.mix(MEGAHIT.out.versions.first())
        }



        GUNZIP_ASSEMBLIES(ch_assembled_contigs)
        ch_versions = ch_versions.mix(GUNZIP_ASSEMBLIES.out.versions)

        ch_assemblies = GUNZIP_ASSEMBLIES.out.gunzip
    }
    else {
        ch_assemblies_split = ch_input_assemblies.branch { _meta, assembly ->
            gzipped: assembly.getExtension() == "gz"
            ungzip: true
        }

        GUNZIP_ASSEMBLYINPUT(ch_assemblies_split.gzipped)
        ch_versions = ch_versions.mix(GUNZIP_ASSEMBLYINPUT.out.versions)

        ch_assemblies = Channel.empty()
        ch_assemblies = ch_assemblies.mix(ch_assemblies_split.ungzip, GUNZIP_ASSEMBLYINPUT.out.gunzip)
    }

    if (!params.skip_quast) {
        QUAST(ch_assemblies)
        ch_versions = ch_versions.mix(QUAST.out.versions.first())
    }

    /*
    ================================================================================
                                    Predict proteins
    ================================================================================
    */

    if (!params.skip_prodigal) {
        PRODIGAL(
            ch_assemblies,
            'gff',
        )
        ch_versions = ch_versions.mix(PRODIGAL.out.versions.first())
    }

    /*
    ================================================================================
                                    Virus identification
    ================================================================================
    */

    if (params.run_virus_identification) {
        VIRUS_IDENTIFICATION(ch_assemblies, ch_genomad_db)
        ch_versions = ch_versions.mix(VIRUS_IDENTIFICATION.out.versions.first())
    }

    /*
    ================================================================================
                                Binning preparation
    ================================================================================
    */

    ch_bin_qc_summary = Channel.empty()

    if (!params.skip_binning || params.ancient_dna) {
        BINNING_PREPARATION(
            ch_assemblies,
            ch_short_reads,
        )
        ch_versions = ch_versions.mix(BINNING_PREPARATION.out.bowtie2_version.first())
    }

    /*
    ================================================================================
                                    Ancient DNA
    ================================================================================
    */

    if (params.ancient_dna) {
        ANCIENT_DNA_ASSEMBLY_VALIDATION(BINNING_PREPARATION.out.grouped_mappings)
        ch_versions = ch_versions.mix(ANCIENT_DNA_ASSEMBLY_VALIDATION.out.versions.first())
    }

    /*
    ================================================================================
                                    Binning
    ================================================================================
    */

    if (!params.skip_binning) {

        // Make sure if running aDNA subworkflow to use the damage-corrected contigs for higher accuracy
        if (params.ancient_dna && !params.skip_ancient_damagecorrection) {
            BINNING(
                BINNING_PREPARATION.out.grouped_mappings.join(ANCIENT_DNA_ASSEMBLY_VALIDATION.out.contigs_recalled).map { it -> [it[0], it[4], it[2], it[3]] }
            )
        }
        else {
            BINNING(
                BINNING_PREPARATION.out.grouped_mappings
            )
        }
        ch_versions = ch_versions.mix(BINNING.out.versions)

        if (params.bin_domain_classification) {

            // Make sure if running aDNA subworkflow to use the damage-corrected contigs for higher accuracy
            if (params.ancient_dna && !params.skip_ancient_damagecorrection) {
                ch_assemblies_for_domainclassification = ANCIENT_DNA_ASSEMBLY_VALIDATION.out.contigs_recalled
            }
            else {
                ch_assemblies_for_domainclassification = ch_assemblies
            }

            DOMAIN_CLASSIFICATION(ch_assemblies_for_domainclassification, BINNING.out.bins, BINNING.out.unbinned)
            ch_binning_results_bins = DOMAIN_CLASSIFICATION.out.classified_bins
            ch_binning_results_unbins = DOMAIN_CLASSIFICATION.out.classified_unbins
            ch_versions = ch_versions.mix(DOMAIN_CLASSIFICATION.out.versions)
        }
        else {
            ch_binning_results_bins = BINNING.out.bins.map { meta, bins ->
                def meta_new = meta + [domain: 'unclassified']
                [meta_new, bins]
            }
            ch_binning_results_unbins = BINNING.out.unbinned.map { meta, bins ->
                def meta_new = meta + [domain: 'unclassified']
                [meta_new, bins]
            }
        }

        /*
        * DAS Tool: binning refinement
        */

        ch_binning_results_bins = ch_binning_results_bins.map { meta, bins ->
            def meta_new = meta + [refinement: 'unrefined']
            [meta_new, bins]
        }

        ch_binning_results_unbins = ch_binning_results_unbins.map { meta, bins ->
            def meta_new = meta + [refinement: 'unrefined_unbinned']
            [meta_new, bins]
        }

        // If any two of the binners are both skipped at once, do not run because DAS_Tool needs at least one
        if (params.refine_bins_dastool) {
            ch_prokarya_bins_dastool = ch_binning_results_bins.filter { meta, _bins ->
                meta.domain != "eukarya"
            }

            if (params.ancient_dna) {
                ch_contigs_for_binrefinement = ANCIENT_DNA_ASSEMBLY_VALIDATION.out.contigs_recalled
            }
            else {
                ch_contigs_for_binrefinement = BINNING_PREPARATION.out.grouped_mappings.map { meta, contigs, _bam, _bai -> [meta, contigs] }
            }

            BINNING_REFINEMENT(ch_contigs_for_binrefinement, ch_prokarya_bins_dastool)

            ch_refined_bins = BINNING_REFINEMENT.out.refined_bins
            ch_refined_unbins = BINNING_REFINEMENT.out.refined_unbins
            ch_versions = ch_versions.mix(BINNING_REFINEMENT.out.versions)

            if (params.postbinning_input == 'raw_bins_only') {
                ch_input_for_postbinning_bins = ch_binning_results_bins
                ch_input_for_postbinning_bins_unbins = ch_binning_results_bins.mix(ch_binning_results_unbins)
            }
            else if (params.postbinning_input == 'refined_bins_only') {
                ch_input_for_postbinning_bins = ch_refined_bins
                ch_input_for_postbinning_bins_unbins = ch_refined_bins.mix(ch_refined_unbins)
            }
            else if (params.postbinning_input == 'both') {
                ch_all_bins = ch_binning_results_bins.mix(ch_refined_bins)
                ch_input_for_postbinning_bins = ch_all_bins
                ch_input_for_postbinning_bins_unbins = ch_all_bins.mix(ch_binning_results_unbins).mix(ch_refined_unbins)
            }
        }
        else {
            ch_input_for_postbinning_bins = ch_binning_results_bins
            ch_input_for_postbinning_bins_unbins = ch_binning_results_bins.mix(ch_binning_results_unbins)
        }

        ch_input_for_postbinning = params.exclude_unbins_from_postbinning
            ? ch_input_for_postbinning_bins
            : ch_input_for_postbinning_bins_unbins

        DEPTHS(ch_input_for_postbinning, BINNING.out.metabat2depths, ch_short_reads)
        ch_input_for_binsummary = DEPTHS.out.depths_summary
        ch_versions = ch_versions.mix(DEPTHS.out.versions)

        /*
        * Bin QC subworkflows: for checking bin completeness with either BUSCO, CHECKM, CHECKM2, and/or GUNC
        */

        if (!params.skip_binqc) {
            BIN_QC(ch_input_for_postbinning)

            ch_bin_qc_summary = BIN_QC.out.qc_summary
            ch_versions = ch_versions.mix(BIN_QC.out.versions)
        }

        ch_quast_bins_summary = Channel.empty()
        if (!params.skip_quast) {
            ch_input_for_quast_bins = ch_input_for_postbinning
                .groupTuple()
                .map { meta, bins ->
                    def new_bins = bins.flatten()
                    [meta, new_bins]
                }

            QUAST_BINS(ch_input_for_quast_bins)
            ch_versions = ch_versions.mix(QUAST_BINS.out.versions.first())
            ch_quast_bin_summary = QUAST_BINS.out.quast_bin_summaries.collectFile(keepHeader: true) { meta, summary ->
                ["${meta.id}.tsv", summary]
            }
            QUAST_BINS_SUMMARY(ch_quast_bin_summary.collect())
            ch_quast_bins_summary = QUAST_BINS_SUMMARY.out.summary
        }

        /*
         * CAT: Bin Annotation Tool (BAT) are pipelines for the taxonomic classification of long DNA sequences and metagenome assembled genomes (MAGs/bins)
         */
        ch_cat_db = Channel.empty()
        if (params.cat_db) {
            CAT_DB(ch_cat_db_file)
            ch_cat_db = CAT_DB.out.db
        }
        else if (params.cat_db_generate) {
            CAT_DB_GENERATE()
            ch_cat_db = CAT_DB_GENERATE.out.db
        }
        CAT(
            ch_input_for_postbinning,
            ch_cat_db,
        )
        // Group all classification results for each sample in a single file
        ch_cat_summary = CAT.out.tax_classification_names.collectFile(keepHeader: true) { meta, classification ->
            ["${meta.id}.txt", classification]
        }
        // Group all classification results for the whole run in a single file
        CAT_SUMMARY(
            ch_cat_summary.collect()
        )
        ch_versions = ch_versions.mix(CAT.out.versions.first())
        ch_versions = ch_versions.mix(CAT_SUMMARY.out.versions)

        // If CAT is not run, then the CAT global summary should be an empty channel
        if (params.cat_db_generate || params.cat_db) {
            ch_cat_global_summary = CAT_SUMMARY.out.combined
        }
        else {
            ch_cat_global_summary = Channel.empty()
        }

        /*
         * GTDB-tk: taxonomic classifications using GTDB reference
         */

        if (!params.skip_gtdbtk) {

            ch_gtdbtk_summary = Channel.empty()
            if (gtdb) {

                ch_gtdb_bins = ch_input_for_postbinning.filter { meta, _bins ->
                    meta.domain != "eukarya"
                }

                GTDBTK(
                    ch_gtdb_bins,
                    ch_bin_qc_summary,
                    gtdb,
                    gtdb_mash,
                )
                ch_versions = ch_versions.mix(GTDBTK.out.versions.first())
                ch_gtdbtk_summary = GTDBTK.out.summary
            }
        }
        else {
            ch_gtdbtk_summary = Channel.empty()
        }

        if ((!params.skip_binqc) || !params.skip_quast || !params.skip_gtdbtk) {
            BIN_SUMMARY(
                ch_input_for_binsummary,
                ch_bin_qc_summary.ifEmpty([]),
                ch_quast_bins_summary.ifEmpty([]),
                ch_gtdbtk_summary.ifEmpty([]),
                ch_cat_global_summary.ifEmpty([]),
                params.binqc_tool,
            )
        }

        /*
         * Prokka: Genome annotation
         */

        if (!params.skip_prokka) {
            ch_bins_for_prokka = ch_input_for_postbinning
                .transpose()
                .map { meta, bin ->
                    def meta_new = meta + [id: bin.getBaseName()]
                    [meta_new, bin]
                }
                .filter { meta, _bin ->
                    meta.domain != "eukarya"
                }

            PROKKA(
                ch_bins_for_prokka,
                [],
                [],
            )
            ch_versions = ch_versions.mix(PROKKA.out.versions.first())
        }

        if (!params.skip_metaeuk && (params.metaeuk_db || params.metaeuk_mmseqs_db)) {
            ch_bins_for_metaeuk = ch_input_for_postbinning
                .transpose()
                .filter { meta, _bin ->
                    meta.domain in ["eukarya", "unclassified"]
                }
                .map { meta, bin ->
                    def meta_new = meta + [id: bin.getBaseName()]
                    [meta_new, bin]
                }

            METAEUK_EASYPREDICT(ch_bins_for_metaeuk, ch_metaeuk_db)
            ch_versions = ch_versions.mix(METAEUK_EASYPREDICT.out.versions)
        }
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'mag_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.fromPath("${workflow.projectDir}/docs/images/mag_logo_mascot_light.png", checkIfExists: true)

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    ch_multiqc_files = ch_multiqc_files.mix(CENTRIFUGE_KREPORT.out.kreport.collect { it[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2.out.report.collect { it[1] }.ifEmpty([]))

    if (!params.skip_quast) {
        ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.report.collect().ifEmpty([]))

        if (!params.skip_binning) {
            ch_multiqc_files = ch_multiqc_files.mix(QUAST_BINS.out.dir.collect().ifEmpty([]))
        }
    }

    if (!params.skip_binning || params.ancient_dna) {
        ch_multiqc_files = ch_multiqc_files.mix(BINNING_PREPARATION.out.bowtie2_assembly_multiqc.collect().ifEmpty([]))
    }

    if (!params.skip_binning && !params.skip_prokka) {
        ch_multiqc_files = ch_multiqc_files.mix(PROKKA.out.txt.collect { it[1] }.ifEmpty([]))
    }

    if (!params.skip_binning && !params.skip_binqc && params.binqc_tool == 'busco') {
        ch_multiqc_files = ch_multiqc_files.mix(BIN_QC.out.multiqc_files.collect().ifEmpty([]))
    }


    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
