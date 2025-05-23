/*
 * CAT/BAT/RAT: tools for taxonomic classification of contigs and metagenome-assembled genomes (MAGs) and for taxonomic profiling of metagenomes
 */
include { CATPACK_ADDNAMES as CATPACK_ADDNAMES_BINS   } from '../../modules/nf-core/catpack/addnames/main'
include { CATPACK_ADDNAMES as CATPACK_ADDNAMES_UNBINS } from '../../modules/nf-core/catpack/addnames/main'
include { CATPACK_BINS                                } from '../../modules/nf-core/catpack/bins/main'
include { CATPACK_CONTIGS                             } from '../../modules/nf-core/catpack/contigs/main'
include { CATPACK_DOWNLOAD                            } from '../../modules/nf-core/catpack/download/main'
include { CATPACK_PREPARE                             } from '../../modules/nf-core/catpack/prepare/main'
include { CATPACK_SUMMARISE                           } from '../../modules/nf-core/catpack/summarise/main'
include { CSVTK_CONCAT as CONCAT_CAT_BINS             } from '../../modules/nf-core/csvtk/concat/main'
include { UNTAR as CAT_DB_UNTAR                       } from '../../modules/nf-core/untar/main'


workflow CATPACK {
    take:
    ch_bins   // channel: [ val(meta), [bins] ]
    ch_unbins // channel: [ val(meta), [unbins] ]

    main:
    ch_versions = Channel.empty()

    /*
    ========================================
     * Database setup
    ========================================
     */

    if (params.cat_db) {
        if (params.cat_db.endsWith('.tar.gz')) {
            CAT_DB_UNTAR([[id: 'cat_db'], file(params.cat_db, checkIfExists: true)])
            ch_cat_db_dir = CAT_DB_UNTAR.out.untar
        }
        else {
            ch_cat_db_dir = Channel.fromPath(params.cat_db, checkIfExists: true, type: 'dir')
                .map { dir -> [[id: 'cat_db'], dir] }
        }

        ch_cat_db = ch_cat_db_dir.multiMap { meta, dir ->
            db: [meta, file(dir / 'db', checkIfExists: true)]
            taxonomy: [meta, file(dir / 'tax', checkIfExists: true)]
        }
    }
    else {
        // download and build the database
        CATPACK_DOWNLOAD([[id: 'cat_db_nr'], 'nr'])
        CATPACK_PREPARE(
            CATPACK_DOWNLOAD.out.fasta,
            CATPACK_DOWNLOAD.out.names.map { _meta, names -> names },
            CATPACK_DOWNLOAD.out.nodes.map { _meta, nodes -> nodes },
            CATPACK_DOWNLOAD.out.acc2tax.map { _meta, acc2tax -> acc2tax },
        )

        ch_cat_db = CATPACK_PREPARE.out
    }

    /*
    =========================================
     * Bin taxonomic classification
    =========================================
     */

    CATPACK_BINS(ch_bins, ch_cat_db.db, ch_cat_db.taxonomy, [[:], []], [[:], []], '.fa')

    ch_bin_classification = CATPACK_BINS.out.bin2classification
        .map { _meta, summary -> [[id: 'bat_bins'], summary] }
        .groupTuple()

    CONCAT_CAT_BINS(ch_bin_classification, 'tsv', 'tsv')

    CATPACK_ADDNAMES_BINS(CONCAT_CAT_BINS.out.csv, ch_cat_db.taxonomy)

    ch_versions = ch_versions.mix(
        CATPACK_BINS.out.versions.first(),
        CONCAT_CAT_BINS.out.versions,
        CATPACK_ADDNAMES_BINS.out.versions,
    )

    /*
    =========================================
     * Unbinned data taxonomic classification
    =========================================
     */

    CATPACK_CONTIGS(ch_unbins, ch_cat_db.db, ch_cat_db.taxonomy, [[:], []], [[:], []])

    CATPACK_ADDNAMES_UNBINS(CATPACK_CONTIGS.out.contig2classification, ch_cat_db.taxonomy)

    ch_unbin_classification = CATPACK_ADDNAMES_UNBINS.out.txt
        .join(ch_unbins)
        .multiMap { meta, names, contigs ->
            names: [meta, names]
            contigs: [meta, contigs]
        }

    CATPACK_SUMMARISE(ch_unbin_classification.names, ch_unbin_classification.contigs)

    ch_versions = ch_versions.mix(
        CATPACK_CONTIGS.out.versions.first(),
        CATPACK_ADDNAMES_UNBINS.out.versions.first(),
        CATPACK_SUMMARISE.out.versions.first(),
    )

    emit:
    summary  = CATPACK_ADDNAMES_BINS.out.txt.map { _meta, summary -> summary }
    versions = ch_versions
}
