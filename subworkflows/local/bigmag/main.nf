/*
 * BUSCO/CheckM2: Alternative workflow to run the missing QC tool to create the file for BIgMAG
 */

include { BUSCO_BUSCO                      } from '../../../modules/nf-core/busco/busco/main'
include { CHECKM2_DATABASEDOWNLOAD         } from '../../../modules/nf-core/checkm2/databasedownload/main'
include { CHECKM2_PREDICT                  } from '../../../modules/nf-core/checkm2/predict/main'
include { UNTAR as BUSCO_UNTAR             } from '../../../modules/nf-core/untar/main'
include { BIGMAG_SUMMARY                   } from '../../../modules/local/bigmag_summary/main'
include { GUNC_DOWNLOADDB                  } from '../../../modules/nf-core/gunc/downloaddb/main'
include { GUNC_RUN                         } from '../../../modules/nf-core/gunc/run/main'
include { CONCAT_BIGMAG                    } from '../../../modules/local/concat_bigmag/main'

workflow BIGMAG {
    take:
    ch_bins // [ [ meta] , fasta ], input bins (mandatory)
    summary // Output from BIN_SUMMARY

    main:
    ch_input_bins_for_qc = ch_bins.transpose()
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    /*
    ================================
     * Setup databases
    ================================
     */

    if (params.busco_db) {
        ch_busco_db = file(params.busco_db, checkIfExists: true)
    }
    else {
        ch_busco_db = []
    }

    if (params.checkm2_db) {
        ch_checkm2_db = [[:], file(params.checkm2_db, checkIfExists: true)]
    }
    else if (params.binqc_tool == 'busco') {
        CHECKM2_DATABASEDOWNLOAD(params.checkm2_db_version)
        ch_checkm2_db = CHECKM2_DATABASEDOWNLOAD.out.database
    }
    else {
        ch_checkm2_db = []
    }

    if (params.gunc_db) {
        ch_gunc_db = file(params.gunc_db, checkIfExists: true)
    }
    else {
        ch_gunc_db = Channel.empty()
    }

    /*
    ================================
     * Run QC tools
    ================================
     */

    if (params.binqc_tool == "checkm" || params.binqc_tool == "checkm2") {
        /*
         * BUSCO
         */

        // Prepare database object depending on type
        // If directory, assumes sub structure of `params.busco_db/lineages/<taxa>_odb(10|12)`
        if (ch_busco_db && ch_busco_db.extension in ['gz', 'tgz']) {
            BUSCO_UNTAR([[id: ch_busco_db.getSimpleName()], ch_busco_db])
            ch_busco_db = BUSCO_UNTAR.out.untar.map { it[1] }
        }
        else if (ch_busco_db && ch_busco_db.isDirectory()) {
            ch_busco_db = ch_busco_db
        }
        else {
            ch_busco_db = []
        }

        BUSCO_BUSCO(ch_bins, 'genome', params.busco_db_lineage, ch_busco_db, [], params.busco_clean)

        ch_alt_summary = BUSCO_BUSCO.out.batch_summary
            .map { _meta, summary -> [[id: 'busco'], summary] }
            .groupTuple()
        ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(
            BUSCO_BUSCO.out.short_summaries_txt.map { it[1] }.flatten()
        )
    }
    else if (params.binqc_tool == "busco") {
        /*
         * CheckM2
         */
        CHECKM2_PREDICT(ch_input_bins_for_qc.groupTuple(), ch_checkm2_db)

        ch_alt_summary = CHECKM2_PREDICT.out.checkm2_tsv
            .map { _meta, summary -> [[id: 'checkm2'], summary] }
            .groupTuple()
        ch_versions = ch_versions.mix(CHECKM2_PREDICT.out.versions.first())
    }

 if (params.run_gunc) {
    ch_gunc_summary = Channel.fromPath( "${params.outdir}/GenomeBinning/QC/gunc_summary.tsv" )
    }

 else if (!params.run_gunc) {
        /*
         * GUNC
         */
        ch_input_bins_for_gunc = ch_bins.filter { meta, _bins ->
            meta.domain != "eukarya"
        }

        if (params.gunc_db) {
            ch_db_for_gunc = ch_gunc_db
        }
        else {
            ch_db_for_gunc = GUNC_DOWNLOADDB(params.gunc_database_type).db
            ch_versions.mix(GUNC_DOWNLOADDB.out.versions)
        }

        GUNC_RUN(ch_input_bins_for_gunc, ch_db_for_gunc)
        ch_versions.mix(GUNC_RUN.out.versions)

        // Make sure to keep directory in sync with modules.conf
        ch_gunc_summary = GUNC_RUN.out.maxcss_level_tsv
            .map { _meta, gunc_summary -> gunc_summary }
            .collectFile(
                name: "gunc_summary.tsv",
                keepHeader: true,
                storeDir: "${params.outdir}/GenomeBinning/QC/",
            )
    }

    // Combine summaries
    CONCAT_BIGMAG (ch_alt_summary, 'tsv', 'tsv')
    ch_missing_summary = CONCAT_BIGMAG.out.csv.map { _meta, summary -> summary }
    ch_versions = ch_versions.mix(CONCAT_BIGMAG.out.versions)

    BIGMAG_SUMMARY(summary, ch_gunc_summary, ch_missing_summary, params.binqc_tool)
    ch_versions = ch_versions.mix(BIGMAG_SUMMARY.out.versions)

    emit:
    bigmag_summary  = BIGMAG_SUMMARY.out.bigmag_summary
    multiqc_files   = ch_multiqc_files
    versions        = ch_versions
}
