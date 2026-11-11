include { METABINNER_KMER       } from '../../../modules/local/metabinner_kmer/main.nf'
include { METABINNER_METABINNER } from '../../../modules/local/metabinner_metabinner/main.nf'
include { METABINNER_BINS       } from '../../../modules/local/metabinner_bins/main.nf'

workflow BINNING_METABINNER {
    take:
    ch_input // [val(meta), path(fasta), path(depth)] (mandatory)

    main:
    ch_versions = channel.empty()

    ch_assembly = ch_input.map { meta, assembly, _depths -> [meta, assembly] }

    // produce k-mer composition table
    METABINNER_KMER(ch_assembly)
    ch_versions = ch_versions.mix(METABINNER_KMER.out.versions)

    // binning
    ch_metabinner_input = ch_input
        .join(METABINNER_KMER.out.composition_profile)
        .map { meta, assembly, depths, kmer -> [meta, assembly, kmer, depths] }
    METABINNER_METABINNER(ch_metabinner_input)
    ch_versions = ch_versions.mix(METABINNER_METABINNER.out.versions)

    // extract bin sequences
    METABINNER_BINS(ch_assembly.join(METABINNER_METABINNER.out.membership))
    ch_versions = ch_versions.mix(METABINNER_BINS.out.versions)

    emit:
    unbinned = METABINNER_BINS.out.unbinned
    bins     = METABINNER_BINS.out.bins
    versions = ch_versions // channel: [ versions.yml ]
}
