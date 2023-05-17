include { CAT_DB          } from '../../modules/local/cat_db'
include { CAT_DB_GENERATE } from '../../modules/local/cat_db_generate'
include { CAT             } from '../../modules/local/cat'
include { CAT_SUMMARY     } from "../../modules/local/cat_summary"
include { GTDBTK          } from '../../subworkflows/local/gtdbtk'

workflow BIN_TAXONOMY {
    take:
    bins_unbins
    busco_summary
    checkm_summary

    main:
    ch_versions = Channel.empty()

    if(params.cat_db){
        ch_cat_db_file = Channel
            .value(file( "${params.cat_db}" ))
    } else {
        ch_cat_db_file = Channel.empty()
    }

    gtdb = params.skip_binqc ? false : params.gtdb
    if (gtdb) {
        ch_gtdb = Channel
            .value(file( "${gtdb}" ))
    } else {
        ch_gtdb = Channel.empty()
    }

    ch_cat_db = Channel.empty()
    if (params.cat_db){
        CAT_DB ( ch_cat_db_file )
        ch_cat_db = CAT_DB.out.db
    } else if (params.cat_db_generate){
        CAT_DB_GENERATE ()
        ch_cat_db = CAT_DB_GENERATE.out.db
    }

    bins_unbins_transposed = bins_unbins.transpose()

    /*
    * CAT: Bin Annotation Tool (BAT) are pipelines for the taxonomic classification of long DNA sequences and metagenome assembled genomes (MAGs/bins)
    */

    CAT (
        bins_unbins_transposed,
        ch_cat_db
    )
    CAT_SUMMARY(
        CAT.out.tax_classification.collect()
    )
    ch_versions = ch_versions.mix(CAT.out.versions.first())
    ch_versions = ch_versions.mix(CAT_SUMMARY.out.versions)

    /*
    * GTDB-tk: taxonomic classifications using GTDB reference
    */
    ch_gtdbtk_summary = Channel.empty()
    if ( gtdb ){
        GTDBTK (
            bins_unbins_transposed,
            ch_busco_summary,
            ch_checkm_summary,
            ch_gtdb
        )
        ch_versions = ch_versions.mix(GTDBTK.out.versions.first())
        ch_gtdbtk_summary = GTDBTK.out.summary
    }

    emit:
    gtdbtk_summary = ch_gtdbtk_summary
    versions       = ch_versions
}
