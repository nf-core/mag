/*
 * BUSCO/CheckM/CheckM2/GUNC: Quantitative measures for the assessment of genome assembly
 */

include { ARIA2 as ARIA2_UNTAR             } from '../../modules/nf-core/aria2/main'
include { BUSCO_BUSCO                      } from '../../modules/nf-core/busco/busco/main'
include { CHECKM2_DATABASEDOWNLOAD         } from '../../modules/nf-core/checkm2/databasedownload/main'
include { CHECKM_QA                        } from '../../modules/nf-core/checkm/qa/main'
include { CHECKM_LINEAGEWF                 } from '../../modules/nf-core/checkm/lineagewf/main'
include { CHECKM2_PREDICT                  } from '../../modules/nf-core/checkm2/predict/main'
include { COMBINE_TSV as COMBINE_BINQC_TSV } from '../../modules/local/combine_tsv'
include { GUNC_DOWNLOADDB                  } from '../../modules/nf-core/gunc/downloaddb/main'
include { GUNC_RUN                         } from '../../modules/nf-core/gunc/run/main'
include { GUNC_MERGECHECKM                 } from '../../modules/nf-core/gunc/mergecheckm/main'
include { UNTAR as BUSCO_UNTAR             } from '../../modules/nf-core/untar'


workflow BIN_QC {
    take:
    ch_bins // [ [ meta] , fasta ], input bins (mandatory)

    main:
    qc_summary = []
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
        busco_lineage = params.busco_auto_lineage_prok ? 'auto_prok' : 'auto'

        if (!ch_busco_db.isEmpty()) {
            if (ch_busco_db.extension in ['gz', 'tgz']) {
                BUSCO_UNTAR([[id: 'busco_db'], ch_busco_db])

                if (ch_busco_db.getSimpleName().contains('odb')) {
                    busco_lineage = ch_busco_db.getSimpleName()
                }
                ch_busco_db = BUSCO_UNTAR.out.untar.map { it[1] }
            }
            else if (ch_busco_db.isDirectory()) {
                if (ch_busco_db.name.matches(/odb\d+$/)) {
                    busco_lineage = ch_busco_db.name
                }
            }
        }

        BUSCO_BUSCO(ch_bins, 'genome', busco_lineage, ch_busco_db, [])

        COMBINE_BINQC_TSV(BUSCO_BUSCO.out.batch_summary.map { it[1] }.collect())

        qc_summary = COMBINE_BINQC_TSV.out.combined
        ch_versions = ch_versions.mix(
            BUSCO_BUSCO.out.versions.first(),
            COMBINE_BINQC_TSV.out.versions
        )
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

        COMBINE_BINQC_TSV(CHECKM_QA.out.output.collect { summary -> summary[1] })

        qc_summary = COMBINE_BINQC_TSV.out.combined
        ch_versions = ch_versions.mix(
            CHECKM_QA.out.versions.first(),
            COMBINE_BINQC_TSV.out.versions
        )
    }
    else if (params.binqc_tool == "checkm2") {
        /*
         * CheckM2
         */
        CHECKM2_PREDICT(ch_input_bins_for_qc.groupTuple(), ch_checkm2_db)

        COMBINE_BINQC_TSV(CHECKM2_PREDICT.out.checkm2_tsv.collect { summary -> summary[1] })

        qc_summary = COMBINE_BINQC_TSV.out.combined
        ch_versions = ch_versions.mix(
            CHECKM2_PREDICT.out.versions.first(),
            COMBINE_BINQC_TSV.out.versions
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
        GUNC_RUN.out.maxcss_level_tsv
            .map { _meta, gunc_summary -> gunc_summary }
            .collectFile(
                name: "gunc_summary.tsv",
                keepHeader: true,
                storeDir: "${params.outdir}/GenomeBinning/QC/"
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
                    storeDir: "${params.outdir}/GenomeBinning/QC/"
                )
        }
    }

    emit:
    qc_summary    = qc_summary
    multiqc_files = ch_multiqc_files
    versions      = ch_versions
}
