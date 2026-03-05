/*
 * CAT/BAT/RAT: tools for taxonomic classification of contigs and metagenome-assembled genomes (MAGs) and for taxonomic profiling of metagenomes
 */
include { CATPACK_ADDNAMES as CATPACK_ADDNAMES_BINS     } from '../../../modules/nf-core/catpack/addnames/main'
include { CATPACK_ADDNAMES as CATPACK_ADDNAMES_UNBINS   } from '../../../modules/nf-core/catpack/addnames/main'
include { CATPACK_BINS                                  } from '../../../modules/nf-core/catpack/bins/main'
include { CATPACK_CONTIGS as CATPACK_UNBINS             } from '../../../modules/nf-core/catpack/contigs/main'
include { CATPACK_DOWNLOAD                              } from '../../../modules/nf-core/catpack/download/main'
include { CATPACK_PREPARE                               } from '../../../modules/nf-core/catpack/prepare/main'
include { CATPACK_SUMMARISE as CATPACK_SUMMARISE_BINS   } from '../../../modules/nf-core/catpack/summarise/main'
include { CATPACK_SUMMARISE as CATPACK_SUMMARISE_UNBINS } from '../../../modules/nf-core/catpack/summarise/main'
include { UNTAR as CAT_DB_UNTAR                         } from '../../../modules/nf-core/untar/main'


workflow CATPACK {
    take:
    ch_bins   // [val(meta), path(fasta)]
    ch_unbins // [val(meta), path(fasta)]

    main:
    ch_versions = channel.empty()

    /*
    ========================================
     * Database setup
    ========================================
     */

    if (params.cat_db) {
        if (params.cat_db.endsWith('.tar.gz')) {
            CAT_DB_UNTAR([[id: 'cat_db'], file(params.cat_db, checkIfExists: true)])
            ch_versions = ch_versions.mix(CAT_DB_UNTAR.out.versions)

            ch_cat_db_dir = CAT_DB_UNTAR.out.untar
        }
        else {
            ch_cat_db_dir = channel.fromPath(params.cat_db, checkIfExists: true, type: 'dir')
                .map { dir -> [[id: 'cat_db'], dir] }
                .first()
        }

        ch_cat_db = ch_cat_db_dir.multiMap { meta, dir ->
            db: [meta, file(dir / 'db', checkIfExists: true)]
            taxonomy: [meta, file(dir / 'tax', checkIfExists: true)]
        }
    }
    else {
        // download and build the database
        log.warn("[nf-core/mag]: Downloading CAT-nr database - this is very large and take a long time!")
        CATPACK_DOWNLOAD([[id: 'cat_db_nr'], 'nr'])
        ch_versions = ch_versions.mix(CATPACK_DOWNLOAD.out.versions)

        CATPACK_PREPARE(
            CATPACK_DOWNLOAD.out.fasta,
            CATPACK_DOWNLOAD.out.names.map { _meta, names -> names },
            CATPACK_DOWNLOAD.out.nodes.map { _meta, nodes -> nodes },
            CATPACK_DOWNLOAD.out.acc2tax.map { _meta, acc2tax -> acc2tax },
        )
        ch_versions = ch_versions.mix(CATPACK_PREPARE.out.versions)

        ch_cat_db = CATPACK_PREPARE.out
    }

    /*
    =========================================
     * Bin taxonomic classification
    =========================================
     */

    CATPACK_BINS(
        ch_bins,
        ch_cat_db.db,
        ch_cat_db.taxonomy,
        [[:], []],
        [[:], []],
        '.fa',
    )
    ch_versions = ch_versions.mix(CATPACK_BINS.out.versions)

    CATPACK_ADDNAMES_BINS(CATPACK_BINS.out.bin2classification, ch_cat_db.taxonomy)
    ch_versions = ch_versions.mix(CATPACK_ADDNAMES_BINS.out.versions)

    bin_summary = CATPACK_ADDNAMES_BINS.out.txt
        .map { _meta, summary -> summary }
        .collectFile(
            name: 'bat_summary.tsv',
            storeDir: "${params.outdir}/Taxonomy/CAT/",
            keepHeader: true,
        )

    if (!params.cat_allow_unofficial_lineages) {
        CATPACK_SUMMARISE_BINS(CATPACK_ADDNAMES_BINS.out.txt, [[:], []])
        ch_versions = ch_versions.mix(CATPACK_SUMMARISE_BINS.out.versions)
    }

    /*
    =========================================
     * Unbinned data taxonomic classification
    =========================================
     */

    if (params.cat_classify_unbinned) {
        CATPACK_UNBINS(
            ch_unbins,
            ch_cat_db.db,
            ch_cat_db.taxonomy,
            [[:], []],
            [[:], []],
        )
        ch_versions = ch_versions.mix(CATPACK_UNBINS.out.versions)

        CATPACK_ADDNAMES_UNBINS(CATPACK_UNBINS.out.contig2classification, ch_cat_db.taxonomy)
        ch_versions = ch_versions.mix(CATPACK_ADDNAMES_UNBINS.out.versions)

        if (!params.cat_allow_unofficial_lineages) {
            ch_unbin_classification = CATPACK_ADDNAMES_UNBINS.out.txt
                .join(ch_unbins)
                .multiMap { meta, names, contigs ->
                    names: [meta, names]
                    contigs: [meta, contigs]
                }

            CATPACK_SUMMARISE_UNBINS(ch_unbin_classification.names, ch_unbin_classification.contigs)
            ch_versions = ch_versions.mix(CATPACK_SUMMARISE_UNBINS.out.versions)
        }
    }

    emit:
    summary  = bin_summary
    versions = ch_versions
}
