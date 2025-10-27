include { METABINNER_KMER           } from '../../../modules/local/metabinner_kmer/main.nf'
include { METABINNER_TOOSHORT       } from '../../../modules/local/metabinner_tooshort/main.nf'
include { METABINNER_METABINNER     } from '../../../modules/local/metabinner_metabinner/main.nf'
include { METABINNER_BINS           } from '../../../modules/local/metabinner_bins/main.nf'

workflow BINNING_METABINNER {

    take:
    ch_input // channel (mandatory): [ val(meta), path(fasta), path(depth) ] (fasta: raw contigs from assembly)

    main:
    ch_versions = Channel.empty()

    // produce k-mer composition table
    METABINNER_KMER(
        ch_input
            .map { meta, assembly, depths ->
                [meta, assembly]
            }
    )
    ch_versions = ch_versions.mix(METABINNER_KMER.out.versions)

    // extract contigs over length threshold
    METABINNER_TOOSHORT(
        ch_input
            .map { meta, assembly, depths ->
                [meta, assembly]
            }
    )
    ch_versions = ch_versions.mix(METABINNER_TOOSHORT.out.versions)

    // binning
    ch_metabinner_input =
        METABINNER_TOOSHORT.out.sizefiltered
        .join(METABINNER_KMER.out.composition_profile)
        .join(ch_input.map { meta, assembly, depths -> [meta, depths] } )
    METABINNER_METABINNER(ch_metabinner_input)
    ch_versions = ch_versions.mix(METABINNER_METABINNER.out.versions)

    // extract bin sequences
    METABINNER_BINS(
        ch_input.map { meta, assembly, depths -> [meta, assembly] }
            .join(METABINNER_METABINNER.out.membership)
    )
    ch_versions = ch_versions.mix(METABINNER_BINS.out.versions)

    emit:
    tooshort            = METABINNER_BINS.out.tooshort
    unbinned            = METABINNER_BINS.out.unbinned
    bins                = METABINNER_BINS.out.bins
    versions            = ch_versions // channel: [ versions.yml ]
}
