/*
 * BUSCO/CheckM/CheckM2/GUNC: Quantitative measures for the assessment of genome assembly
 */

include { BUSCO_DB_PREPARATION              } from '../../modules/local/busco_db_preparation'
include { BUSCO                             } from '../../modules/local/busco'
include { BUSCO_SAVE_DOWNLOAD               } from '../../modules/local/busco_save_download'
include { BUSCO_SUMMARY                     } from '../../modules/local/busco_summary'
include { CHECKM_QA                         } from '../../modules/nf-core/checkm/qa/main'
include { CHECKM_LINEAGEWF                  } from '../../modules/nf-core/checkm/lineagewf/main'
include { CHECKM2_PREDICT                   } from '../../modules/nf-core/checkm2/predict/main'
include { COMBINE_TSV as COMBINE_BINQC_TSV  } from '../../modules/local/combine_tsv'
include { GUNC_DOWNLOADDB                   } from '../../modules/nf-core/gunc/downloaddb/main'
include { GUNC_RUN                          } from '../../modules/nf-core/gunc/run/main'
include { GUNC_MERGECHECKM                  } from '../../modules/nf-core/gunc/mergecheckm/main'


workflow BIN_QC {
    take:
    ch_bins        // [ [ meta] , fasta           ], input bins (mandatory)
    ch_checkm_db   // [ db                        ], presupplied CheckM database (optional)
    ch_checkm2_db  // [ [meta]  , db              ], presupplied CheckM2 database (optional)
    ch_busco_db    // [ [meta]  , db              ], presupplied BUSCO database (optional)
    ch_gunc_db     // [ db                        ], presupplied GUNC database (optional)

    main:
    qc_summary = []
    ch_input_bins_for_qc = ch_bins.transpose()
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


    if (params.binqc_tool == "busco") {
        /*
         * BUSCO
         */
        if (!ch_busco_db.isEmpty()) {
            if (ch_busco_db.extension in ['gz', 'tgz']) {
                // Expects to be tar.gz!
                BUSCO_DB_PREPARATION(ch_busco_db)
                ch_db_for_busco = BUSCO_DB_PREPARATION.out.db.map { meta, db ->
                    [[id: meta, lineage: 'Y'], db]
                }
            }
            else if (ch_busco_db.isDirectory()) {
                // Set meta to match expected channel cardinality for BUSCO
                ch_db_for_busco = Channel
                    .of(ch_busco_db)
                    .map { db ->
                        def basename = db.getBaseName()
                        def lineage = basename.contains('odb10') ? 'Y' : 'N'
                        [[id: basename, lineage: lineage], db]
                    }
                    .collect()
            }
        }
        else {
            // Set BUSCO database to empty to allow for --auto-lineage
            ch_db_for_busco = Channel
                .of([[lineage: ''], []])
                .collect()
        }

        if (params.save_busco_db) {
            // publish files downloaded by Busco
            ch_downloads = BUSCO.out.busco_downloads
                .groupTuple()
                .map { _lin, downloads -> downloads[0] }
                .toSortedList()
                .flatten()
            BUSCO_SAVE_DOWNLOAD(ch_downloads)
        }

        BUSCO(ch_input_bins_for_qc, ch_db_for_busco)

        BUSCO_SUMMARY(
            BUSCO.out.summary_domain.map { it[1] }.collect().ifEmpty([]),
            BUSCO.out.summary_specific.map { it[1] }.collect().ifEmpty([]),
            BUSCO.out.failed_bin.map { it[1] }.collect().ifEmpty([])
        )

        ch_multiqc_files = ch_multiqc_files.mix(
            BUSCO.out.summary_domain.mix(BUSCO.out.summary_specific).map{ it[1] }
        )
        qc_summary = BUSCO_SUMMARY.out.summary
        ch_versions = ch_versions.mix(BUSCO.out.versions.first())
    }
    else if (params.binqc_tool == "checkm") {
        /*
         * CheckM
         */
        ch_bins_for_checkmlineagewf = ch_input_bins_for_qc
            .groupTuple()
            .filter { meta, _bins ->
                    meta.domain != "eukarya"
                }
            .multiMap { meta, fa ->
                reads: [meta, fa]
                ext: fa.extension.unique().join("")  // the pipeline ensures that all bins will have the same extension
            }

        CHECKM_LINEAGEWF(ch_bins_for_checkmlineagewf.reads, ch_bins_for_checkmlineagewf.ext, ch_checkm_db)
        ch_versions = ch_versions.mix(CHECKM_LINEAGEWF.out.versions.first())

        ch_checkmqa_input = CHECKM_LINEAGEWF.out.checkm_output
            .join(CHECKM_LINEAGEWF.out.marker_file)
            .map { meta, dir, marker ->
                [meta, dir, marker, []]
            }

        CHECKM_QA(ch_checkmqa_input, [])

        COMBINE_BINQC_TSV(CHECKM_QA.out.output.map { it[1] }.collect())

        qc_summary = COMBINE_BINQC_TSV.out.combined
        ch_versions = ch_versions.mix(CHECKM_QA.out.versions.first())
    }
    else if (params.binqc_tool == "checkm2") {
        /*
         * CheckM2
         */
        CHECKM2_PREDICT(ch_input_bins_for_qc.groupTuple(), ch_checkm2_db)

        COMBINE_BINQC_TSV(CHECKM2_PREDICT.out.checkm2_tsv.map { it[1] }.collect())

        qc_summary = COMBINE_BINQC_TSV.out.combined
        ch_versions = ch_versions.mix(CHECKM2_PREDICT.out.versions.first())
    }

    if (params.run_gunc) {
        /*
         * GUNC
         */
        ch_input_bins_for_gunc = ch_bins
            .filter { meta, _bins ->
                meta.domain != "eukarya"
            }
            .flatMap { meta, bins ->
                bins.collect { bin -> [meta, bin] }
            }

        if ( params.gunc_db ) {
            ch_db_for_gunc = ch_gunc_db
        }
        else {
            ch_db_for_gunc = GUNC_DOWNLOADDB(params.gunc_database_type).db
            ch_versions.mix(GUNC_DOWNLOADDB.out.versions)
        }

        GUNC_RUN(ch_input_bins_for_gunc, ch_db_for_gunc)
        ch_versions.mix(GUNC_RUN.out.versions)

        // Make sure to keep directory in sync with modules.conf
        GUNC_RUN.out.maxcss_level_tsv
            .map{it[1]}
            .collectFile(name: "gunc_summary.tsv", keepHeader: true, storeDir: "${params.outdir}/GenomeBinning/QC/")

        if ( params.binqc_tool == 'checkm' ) {
            ch_input_to_mergecheckm = GUNC_RUN.out.maxcss_level_tsv.combine(CHECKM_QA.out.output, by: 0)

            GUNC_MERGECHECKM(ch_input_to_mergecheckm)
            ch_versions.mix(GUNC_MERGECHECKM.out.versions)

            // Make sure to keep directory in sync with modules.conf
            GUNC_MERGECHECKM.out.tsv
                .map{it[1]}
                .collectFile(name: "gunc_checkm_summary.tsv", keepHeader: true, storeDir: "${params.outdir}/GenomeBinning/QC/")
        }
    }


    emit:
    qc_summary      = qc_summary
    multiqc_files   = ch_multiqc_files
    versions        = ch_versions
}
