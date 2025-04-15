/*
 * BUSCO/CheckM/CheckM2/GUNC: Quantitative measures for the assessment of genome assembly
 */

include { ARIA2 as ARIA2_UNTAR             } from '../../modules/nf-core/aria2/main'
include { BUSCO_BUSCO                      } from '../../modules/nf-core/busco/busco/main'
include { CHECKM2_DATABASEDOWNLOAD         } from '../../modules/nf-core/checkm2/databasedownload/main'
include { CHECKM_QA                        } from '../../modules/nf-core/checkm/qa/main'
include { CHECKM_LINEAGEWF                 } from '../../modules/nf-core/checkm/lineagewf/main'
include { CHECKM2_PREDICT                  } from '../../modules/nf-core/checkm2/predict/main'
include { CSVTK_CONCAT as CONCAT_BINQC_TSV } from '../../modules/nf-core/csvtk/concat/main'
include { GUNC_DOWNLOADDB                  } from '../../modules/nf-core/gunc/downloaddb/main'
include { GUNC_RUN                         } from '../../modules/nf-core/gunc/run/main'
include { GUNC_MERGECHECKM                 } from '../../modules/nf-core/gunc/mergecheckm/main'
include { UNTAR as BUSCO_UNTAR             } from '../../modules/nf-core/untar/main'


workflow BIN_QC {
    take:
    ch_bins // [ [ meta] , fasta ], input bins (mandatory)

    main:
    ch_qc_summary = []
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

    if (params.checkm_db) {
        ch_checkm_db = file(params.checkm_db, checkIfExists: true)
    }
    else if (params.binqc_tool == 'checkm') {
        ARIA2_UNTAR(params.checkm_download_url)
        ch_checkm_db = ARIA2_UNTAR.out.downloaded_file
    }
    else {
        ch_checkm_db = []
    }

    if (params.checkm2_db) {
        ch_checkm2_db = [[:], file(params.checkm2_db, checkIfExists: true)]
    }
    else if (params.binqc_tool == 'checkm2') {
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

    if (params.binqc_tool == "busco") {
        /*
         * BUSCO
         */

        // Prepare database object depending on type
        // If directory, assumes sub structure of `params.busco_db/lineages/<taxa>_odb(10|12)`
        if (ch_busco_db.extension in ['gz', 'tgz']) {
            BUSCO_UNTAR([[id: ch_busco_db.getSimpleName()], ch_busco_db])
            ch_busco_db = BUSCO_UNTAR.out.untar.map { it[1] }
        }
        else if (ch_busco_db.isDirectory()) {
            ch_busco_db = ch_busco_db
        }
        else {
            ch_busco_db = []
        }

        // Specify a specific lineage otherwise just let BUSCO auto-select
        // Warning: if not all lineages in `--busco_db` it will try to auto-download!
        if (params.busco_db_lineage) {
            busco_lineage = params.busco_db_lineage
        }
        else if (params.busco_auto_lineage_prok) {
            busco_lineage = 'auto_prok'
        }
        else {
            busco_lineage = 'auto'
        }

        BUSCO_BUSCO(ch_bins, 'genome', busco_lineage, ch_busco_db, [], params.busco_clean)

        ch_qc_summaries = BUSCO_BUSCO.out.batch_summary
            .map { _meta, summary -> [[id: 'busco'], summary] }
            .groupTuple()
        ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(
            BUSCO_BUSCO.out.short_summaries_txt.map { it[1] }.flatten()
        )
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
                ext: fa.extension.unique().join("")
            }

        CHECKM_LINEAGEWF(ch_bins_for_checkmlineagewf.reads, ch_bins_for_checkmlineagewf.ext, ch_checkm_db)
        ch_versions = ch_versions.mix(CHECKM_LINEAGEWF.out.versions.first())

        ch_checkmqa_input = CHECKM_LINEAGEWF.out.checkm_output
            .join(CHECKM_LINEAGEWF.out.marker_file)
            .map { meta, dir, marker ->
                [meta, dir, marker, []]
            }

        CHECKM_QA(ch_checkmqa_input, [])

        ch_qc_summaries = CHECKM_QA.out.output
            .map { _meta, summary -> [[id: 'checkm'], summary] }
            .groupTuple()
        ch_versions = ch_versions.mix(CHECKM_QA.out.versions.first())
    }
    else if (params.binqc_tool == "checkm2") {
        /*
         * CheckM2
         */
        CHECKM2_PREDICT(ch_input_bins_for_qc.groupTuple(), ch_checkm2_db)

        ch_qc_summaries = CHECKM2_PREDICT.out.checkm2_tsv
            .map { _meta, summary -> [[id: 'checkm2'], summary] }
            .groupTuple()
        ch_versions = ch_versions.mix(CHECKM2_PREDICT.out.versions.first())
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
        GUNC_RUN.out.maxcss_level_tsv
            .map { _meta, gunc_summary -> gunc_summary }
            .collectFile(
                name: "gunc_summary.tsv",
                keepHeader: true,
                storeDir: "${params.outdir}/GenomeBinning/QC/",
            )

        if (params.binqc_tool == 'checkm') {
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

    // Combine QC summaries (same process for all tools)
    CONCAT_BINQC_TSV(ch_qc_summaries, 'tsv', 'tsv')
    ch_qc_summary = CONCAT_BINQC_TSV.out.csv.map { _meta, summary -> summary }
    ch_versions = ch_versions.mix(CONCAT_BINQC_TSV.out.versions)

    emit:
    qc_summary    = ch_qc_summary
    multiqc_files = ch_multiqc_files
    versions      = ch_versions
}
