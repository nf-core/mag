/*
 * CAT: Bin Annotation Tool (BAT) are pipelines for the taxonomic classification of long DNA sequences and metagenome assembled genomes (MAGs/bins)
 */
include { CATPACK_BINS      } from '../../modules/nf-core/catpack/bins/main'
include { CATPACK_CONTIGS   } from '../../modules/nf-core/catpack/contigs/main'

include { CAT_DB            } from '../../modules/local/cat_db'
include { CAT_DB_GENERATE   } from '../../modules/local/cat_db_generate'
include { CAT as CATOLD     } from '../../modules/local/cat'
include { CAT_SUMMARY       } from '../../modules/local/cat_summary'


workflow CAT {
    take:
    ch_bins   // channel: [ val(meta), [bins] ]
    ch_unbins // channel: [ val(meta), [unbins] ]

    main:
    ch_versions = Channel.empty()

    /*
    ========================================
     * Setup database
    ========================================
     */

    if (params.cat_db) {
        ch_cat_db_file = Channel.value(file("${params.cat_db}"))
    }
    else {
        ch_cat_db_file = Channel.empty()
    }

    ch_cat_db = Channel.empty()
    ch_cat_taxonomy = Channel.empty()
    if (params.cat_db) {
        CAT_DB(ch_cat_db_file)
        ch_cat_db = CAT_DB.out.db
        ch_cat_taxonomy = CAT_DB.out.taxonomy
    }
    else if (params.cat_db_generate) {
        CAT_DB_GENERATE()
        ch_cat_db = CAT_DB_GENERATE.out.db
    }

    /*
    =========================================
     * Bin taxonomic classification
    =========================================
     */

    CATPACK_BINS(
        ch_bins,
        ch_cat_db,
        ch_cat_taxonomy,
        [[:], []],
        [[:], []],
        '.fa',
    )
    ch_versions = ch_versions.mix(CATPACK_BINS.out.versions.first())

    /*
    =========================================
     * Unbinned data taxonomic classification
    =========================================
     */

    CATPACK_CONTIGS(
        ch_unbins,
        ch_cat_db,
        ch_cat_taxonomy,
        [[:], []],
        [[:], []],
    )
    ch_versions = ch_versions.mix(CATPACK_CONTIGS.out.versions.first())

    emit:
    versions = ch_versions
}
