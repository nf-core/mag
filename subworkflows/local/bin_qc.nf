/*
 * BUSCO/CheckM/CheckM2: Quantitative measures for the assessment of genome assembly
 */

include { BUSCO_DB_PREPARATION  } from '../../modules/local/busco_db_preparation'
include { BUSCO                 } from '../../modules/local/busco'
include { BUSCO_SAVE_DOWNLOAD   } from '../../modules/local/busco_save_download'
include { BUSCO_SUMMARY         } from '../../modules/local/busco_summary'
include { CHECKM_QA             } from '../../modules/nf-core/checkm/qa/main'
include { CHECKM_LINEAGEWF      } from '../../modules/nf-core/checkm/lineagewf/main'
include { CHECKM2_PREDICT       } from '../../modules/nf-core/checkm2/predict/main'
include { COMBINE_TSV           } from '../../modules/local/combine_tsv'

workflow BIN_QC {
    take:
    bins       // channel: [ val(meta), path(bin) ]
    checkm_db
    checkm2_db
    busco_db

    main:
    ch_versions = Channel.empty()

    if (params.binqc_tool == "busco") {
        // BUSCO workflow
        if (!busco_db.isEmpty()) {
            if (busco_db.extension in ['gz', 'tgz']) {
                // Expects to be tar.gz!
                BUSCO_DB_PREPARATION(busco_db)
                ch_db_for_busco = BUSCO_DB_PREPARATION.out.db.map { meta, db ->
                    [[id: meta, lineage: 'Y'], db]
                }
            }
            else if (busco_db.isDirectory()) {
                // Set meta to match expected channel cardinality for BUSCO
                ch_db_for_busco = Channel
                    .of(busco_db)
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
                .of([])
                .map { empty_db -> [[lineage: ''], []] }
                .collect()
        }

        if (params.save_busco_db) {
            // publish files downloaded by Busco
            ch_downloads = BUSCO.out.busco_downloads
                .groupTuple()
                .map { lin, downloads -> downloads[0] }
                .toSortedList()
                .flatten()
            BUSCO_SAVE_DOWNLOAD(ch_downloads)
        }

        BUSCO(bins, ch_db_for_busco)

        // busco_summary_domain = BUSCO.out.summary_domain.collect()
        // busco_summary_specific = BUSCO.out.summary_specific.collect()
        // busco_failed_bin = BUSCO.out.failed_bin.collect()

        BUSCO_SUMMARY(
            BUSCO.out.summary_domain.map { it[1] }.collect().ifEmpty([]),
            BUSCO.out.summary_specific.map { it[1] }.collect().ifEmpty([]),
            BUSCO.out.failed_bin.map { it[1] }.collect().ifEmpty([])
        )

        multiqc_reports = BUSCO.out.summary_domain.mix(BUSCO.out.summary_specific).map{ it[1] }
        summary  = BUSCO_SUMMARY.out.summary
        ch_versions = ch_versions.mix(BUSCO.out.versions.first())
    }
    else if (params.binqc_tool == "checkm") {
        // CheckM workflow
        ch_bins_for_checkmlineagewf = bins
            .filter { meta, bin ->
                    meta.domain != "eukarya"
                }
            .multiMap { meta, fa ->
                reads: [meta, fa]
                ext: fa.extension.unique().join("")
            }

        CHECKM_LINEAGEWF(ch_bins_for_checkmlineagewf.reads, ch_bins_for_checkmlineagewf.ext, checkm_db)
        ch_versions = ch_versions.mix(CHECKM_LINEAGEWF.out.versions.first())

        ch_checkmqa_input = CHECKM_LINEAGEWF.out.checkm_output
            .join(CHECKM_LINEAGEWF.out.marker_file)
            .map { meta, dir, marker ->
                [meta, dir, marker, []]
            }

        CHECKM_QA(ch_checkmqa_input, [])

        COMBINE_TSV(CHECKM_QA.out.output.map { it[1] }.collect())

        summary = COMBINE_TSV.out.combined
        ch_versions = ch_versions.mix(CHECKM_QA.out.versions.first())
    }
    else if (params.binqc_tool == "checkm2") {
        // CheckM2 workflow
        CHECKM2_PREDICT(bins, checkm2_db)

        COMBINE_TSV(CHECKM2_PREDICT.out.checkm2_tsv.map { it[1] }.collect())

        summary = COMBINE_TSV.out.combined
        ch_versions = ch_versions.mix(CHECKM2_PREDICT.out.versions.first())
    }

    emit:
    summary     = summary
    checkm_tsv  = params.binqc_tool == "checkm" ? CHECKM_QA.out.output : []
    multiqc     = params.binqc_tool == "busco" ? multiqc_reports : []
    versions    = ch_versions
}
