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
include { FIND_CONCATENATE as CONCAT_UNBINS             } from '../../../modules/nf-core/find/concatenate/main'
include { UNTAR as CATPACK_DB_UNTAR                     } from '../../../modules/nf-core/untar/main'

workflow CATPACK {
    take:
    ch_bins   // [val(meta), path(fasta)]
    ch_unbins // [val(meta), path(fasta)]

    main:
    ch_versions = channel.empty()
    ch_publish_db = channel.empty()
    ch_publish_unbins = channel.empty()

    /*
    ========================================
     * Database setup
    ========================================
     */

    if (params.cat_db) {
        if (params.cat_db.name.endsWith('.tar.gz')) {
            CATPACK_DB_UNTAR([[id: 'cat_db'], file(params.cat_db, checkIfExists: true)])
            ch_versions = ch_versions.mix(CATPACK_DB_UNTAR.out.versions)

            ch_cat_db_dir = CATPACK_DB_UNTAR.out.untar
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

        CATPACK_PREPARE(
            CATPACK_DOWNLOAD.out.fasta,
            CATPACK_DOWNLOAD.out.names.map { _meta, names -> names },
            CATPACK_DOWNLOAD.out.nodes.map { _meta, nodes -> nodes },
            CATPACK_DOWNLOAD.out.acc2tax.map { _meta, acc2tax -> acc2tax },
        )

        ch_cat_db = CATPACK_PREPARE.out
        ch_publish_db = params.save_cat_db ? CATPACK_PREPARE.out.db.mix(CATPACK_PREPARE.out.taxonomy) : channel.empty()
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
        '.fa.gz',
    )

    CATPACK_ADDNAMES_BINS(CATPACK_BINS.out.bin2classification, ch_cat_db.taxonomy)

    ch_publish_bins = CATPACK_BINS.out.orf2lca.mix(
        CATPACK_BINS.out.bin2classification,
        CATPACK_BINS.out.log,
        CATPACK_BINS.out.diamond,
        CATPACK_BINS.out.faa,
        CATPACK_BINS.out.gff,
        CATPACK_ADDNAMES_BINS.out.txt,
    )

    bin_summary = CATPACK_ADDNAMES_BINS.out.txt
        .map { _meta, summary -> summary }
        .collectFile(
            name: 'bat_summary.tsv',
            storeDir: "${params.outdir}/Taxonomy/CAT/",
            keepHeader: true,
            sort: 'deep',
        )

    if (!params.cat_allow_unofficial_lineages) {
        CATPACK_SUMMARISE_BINS(CATPACK_ADDNAMES_BINS.out.txt, [[:], []])
        ch_publish_bins = ch_publish_bins.mix(CATPACK_SUMMARISE_BINS.out.txt)
    }

    /*
    =========================================
     * Unbinned data taxonomic classification
    =========================================
     */

    if (params.cat_classify_unbinned) {
        CONCAT_UNBINS(ch_unbins)

        CATPACK_UNBINS(
            CONCAT_UNBINS.out.file_out,
            ch_cat_db.db,
            ch_cat_db.taxonomy,
            [[:], []],
            [[:], []],
        )

        CATPACK_ADDNAMES_UNBINS(CATPACK_UNBINS.out.contig2classification, ch_cat_db.taxonomy)

        ch_publish_unbins = CATPACK_UNBINS.out.orf2lca.mix(
            CATPACK_UNBINS.out.contig2classification,
            CATPACK_UNBINS.out.log,
            CATPACK_UNBINS.out.diamond,
            CATPACK_UNBINS.out.faa,
            CATPACK_UNBINS.out.gff,
            CATPACK_ADDNAMES_UNBINS.out.txt,
        )

        if (!params.cat_allow_unofficial_lineages) {
            ch_unbin_classification = CATPACK_ADDNAMES_UNBINS.out.txt
                .join(CONCAT_UNBINS.out.file_out)
                .multiMap { meta, names, contigs ->
                    names: [meta, names]
                    contigs: [meta, contigs]
                }

            CATPACK_SUMMARISE_UNBINS(ch_unbin_classification.names, ch_unbin_classification.contigs)
            ch_publish_unbins = ch_publish_unbins.mix(CATPACK_SUMMARISE_UNBINS.out.txt)
        }
    }

    emit:
    summary        = bin_summary
    versions       = ch_versions
    publish_db     = ch_publish_db // [val(meta), path(dir)]
    publish_bins   = ch_publish_bins // [val(meta), path(files)]
    publish_unbins = ch_publish_unbins // [val(meta), path(files)]
}
