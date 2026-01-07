/*
 * BUSCO/CheckM/CheckM2/GUNC: Quantitative measures for the assessment of genome assembly
 */

include { BUSCO_BUSCO                   } from '../../../modules/nf-core/busco/busco/main'
include { CHECKM2_DATABASEDOWNLOAD      } from '../../../modules/nf-core/checkm2/databasedownload/main'
include { CHECKM_QA                     } from '../../../modules/nf-core/checkm/qa/main'
include { CHECKM_LINEAGEWF              } from '../../../modules/nf-core/checkm/lineagewf/main'
include { CHECKM2_PREDICT               } from '../../../modules/nf-core/checkm2/predict/main'
include { QSV_CAT as CONCAT_BUSCO_TSV   } from '../../../modules/nf-core/qsv/cat/main'
include { QSV_CAT as CONCAT_CHECKM_TSV  } from '../../../modules/nf-core/qsv/cat/main'
include { QSV_CAT as CONCAT_CHECKM2_TSV } from '../../../modules/nf-core/qsv/cat/main'
include { GUNC_DOWNLOADDB               } from '../../../modules/nf-core/gunc/downloaddb/main'
include { GUNC_RUN                      } from '../../../modules/nf-core/gunc/run/main'
include { GUNC_MERGECHECKM              } from '../../../modules/nf-core/gunc/mergecheckm/main'
include { UNTAR as BUSCO_UNTAR          } from '../../../modules/nf-core/untar/main'
include { UNTAR as CHECKM_UNTAR         } from '../../../modules/nf-core/untar/main'


workflow BIN_QC {
    take:
    ch_bins // [val(meta), path(fasta)], input bins (mandatory)

    main:
    ch_qc_summaries = channel.empty()
    ch_input_bins_for_qc = ch_bins.transpose()
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_busco_final_summaries = channel.empty()
    ch_checkm_final_summaries = channel.empty()
    ch_checkm2_final_summaries = channel.empty()
    ch_gunc_summary = channel.empty()

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

    if (params.checkm_db) {
        ch_checkm_db = file(params.checkm_db, checkIfExists: true)
    }
    else if (params.run_checkm) {
        ch_checkm_db = [[id: 'checkm_db'], file(params.checkm_download_url, checkIfExists: true)]

        CHECKM_UNTAR(ch_checkm_db)
        ch_versions = ch_versions.mix(CHECKM_UNTAR.out.versions)

        ch_checkm_db = CHECKM_UNTAR.out.untar.map { _meta, db -> db }
    }
    else {
        ch_checkm_db = []
    }

    if (params.checkm2_db) {
        ch_checkm2_db = [[:], file(params.checkm2_db, checkIfExists: true)]
    }
    else if (params.run_checkm2) {
        CHECKM2_DATABASEDOWNLOAD(params.checkm2_db_version)
        ch_versions = ch_versions.mix(CHECKM2_DATABASEDOWNLOAD.out.versions)
        ch_checkm2_db = CHECKM2_DATABASEDOWNLOAD.out.database
    }
    else {
        ch_checkm2_db = []
    }

    if (params.gunc_db) {
        ch_gunc_db = file(params.gunc_db, checkIfExists: true)
    }
    else {
        ch_gunc_db = channel.empty()
    }

    /*
    ================================
     * Run QC tools
    ================================
     */
    if (params.run_busco) {
        /*
         * BUSCO
         */

        // Prepare database object depending on type
        // If directory, assumes sub structure of `params.busco_db/lineages/<taxa>_odb(10|12)`
        if (ch_busco_db && ch_busco_db.extension in ['gz', 'tgz']) {
            BUSCO_UNTAR([[id: ch_busco_db.getSimpleName()], ch_busco_db])
            ch_versions = ch_versions.mix(BUSCO_UNTAR.out.versions)
            ch_busco_db = BUSCO_UNTAR.out.untar.map { _meta, db -> db }
        }
        else if (ch_busco_db && ch_busco_db.isDirectory()) {
            ch_busco_db = ch_busco_db
        }
        else {
            ch_busco_db = []
        }

        BUSCO_BUSCO(ch_bins, 'genome', params.busco_db_lineage, ch_busco_db, [], params.busco_clean)
        ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions)

        ch_busco_summaries = BUSCO_BUSCO.out.batch_summary
            .map { _meta, summary -> [[id: 'busco'], summary] }
            .groupTuple()
        ch_multiqc_files = ch_multiqc_files.mix(
            BUSCO_BUSCO.out.short_summaries_txt.map { _meta, summary -> summary }.flatten()
        )

        CONCAT_BUSCO_TSV(ch_busco_summaries, 'rowskey', 'tsv', true)
        ch_busco_final_summaries = ch_busco_final_summaries.mix(
            CONCAT_BUSCO_TSV.out.csv.map { _meta, csv -> csv }
        )
        ch_qc_summaries = ch_qc_summaries.mix(
            CONCAT_BUSCO_TSV.out.csv.splitCsv(header: true, sep: '\t').map { _meta, summary -> [bin_qc_tool: 'busco'] + summary }
        )
    }
    if (params.run_checkm) {
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
                ext: fa.extension.unique().join("")
            }

        CHECKM_LINEAGEWF(ch_bins_for_checkmlineagewf.reads, ch_bins_for_checkmlineagewf.ext, ch_checkm_db)
        ch_versions = ch_versions.mix(CHECKM_LINEAGEWF.out.versions)

        ch_checkmqa_input = CHECKM_LINEAGEWF.out.checkm_output
            .join(CHECKM_LINEAGEWF.out.marker_file)
            .map { meta, dir, marker ->
                [meta, dir, marker, []]
            }

        CHECKM_QA(ch_checkmqa_input, [])
        ch_versions = ch_versions.mix(CHECKM_QA.out.versions)

        ch_checkm_summaries = CHECKM_QA.out.output
            .map { _meta, summary -> [[id: 'checkm'], summary] }
            .groupTuple()
        ch_multiqc_files = ch_multiqc_files.mix(
            CHECKM_QA.out.output.map { _meta, summary -> summary }.flatten()
        )

        CONCAT_CHECKM_TSV(ch_checkm_summaries, 'rowskey', 'tsv', true)
        ch_checkm_final_summaries = ch_checkm_final_summaries.mix(
            CONCAT_CHECKM_TSV.out.csv.map { _meta, csv -> csv }
        )
        ch_qc_summaries = ch_qc_summaries.mix(
            CONCAT_CHECKM_TSV.out.csv.splitCsv(header: true, sep: '\t').map { _meta, summary -> [bin_qc_tool: 'checkm'] + summary }
        )
    }
    if (params.run_checkm2) {
        /*
         * CheckM2
         */
        CHECKM2_PREDICT(ch_input_bins_for_qc.groupTuple(), ch_checkm2_db)
        ch_versions = ch_versions.mix(CHECKM2_PREDICT.out.versions)

        ch_checkm2_summaries = CHECKM2_PREDICT.out.checkm2_tsv
            .map { _meta, summary -> [[id: 'checkm2'], summary] }
            .groupTuple()
        ch_multiqc_files = ch_multiqc_files.mix(
            CHECKM2_PREDICT.out.checkm2_tsv.map { _meta, summary -> summary }.flatten()
        )

        CONCAT_CHECKM2_TSV(ch_checkm2_summaries, 'rowskey', 'tsv', false)
        ch_checkm2_final_summaries = ch_checkm2_final_summaries.mix(
            CONCAT_CHECKM2_TSV.out.csv.map { _meta, csv -> csv }
        )
        ch_qc_summaries = ch_qc_summaries.mix(
            CONCAT_CHECKM2_TSV.out.csv.splitCsv(header: true, sep: '\t').map { _meta, summary -> [bin_qc_tool: 'checkm2'] + summary }
        )
    }

    if (params.run_gunc) {
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
        if (params.run_checkm) {
            ch_input_to_mergecheckm = GUNC_RUN.out.maxcss_level_tsv.combine(CHECKM_QA.out.output, by: 0)

            GUNC_MERGECHECKM(ch_input_to_mergecheckm)
            ch_versions.mix(GUNC_MERGECHECKM.out.versions)

            // Make sure to keep directory in sync with modules.conf
            GUNC_MERGECHECKM.out.tsv
                .map { _meta, gunc_checkm_summary -> gunc_checkm_summary }
                .collectFile(
                    name: "gunc_checkm_summary.tsv",
                    keepHeader: true,
                    storeDir: "${params.outdir}/GenomeBinning/QC/",
                )
        }
    }

    emit:
    qc_summaries    = ch_qc_summaries
    busco_summary   = ch_busco_final_summaries
    checkm_summary  = ch_checkm_final_summaries
    checkm2_summary = ch_checkm2_final_summaries
    gunc_summary    = ch_gunc_summary
    multiqc_files   = ch_multiqc_files
    versions        = ch_versions
}
