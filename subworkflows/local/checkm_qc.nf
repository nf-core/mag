/*
 * CheckM: Quantitative measures for the assessment of genome assembly
 */

include { CHECKM_QA        } from '../../modules/nf-core/checkm/qa/main'
include { CHECKM_LINEAGEWF } from '../../modules/nf-core/checkm/lineagewf/main'

workflow CHECKM_QC {
    take:
    bins       // channel: [ val(meta), path(bin) ]
    checkm_db

    main:
    // todo initialise checkm_db as optional if user supplied
    CHECKM_LINEAGEWF ( bins, checkm_db )

    // pass output of LINEAGEWF to QA - check snakemake


    emit:
    checkm_tsv = CHECKM_QA.out.output
}
