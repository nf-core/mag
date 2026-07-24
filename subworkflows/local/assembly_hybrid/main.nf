// MODULES
include { SPADES as METASPADESHYBRID } from '../../../modules/nf-core/spades/main'

workflow HYBRID_ASSEMBLY {
    take:
    ch_short_reads_spades // [val(meta), path(fastq1), path(fastq2)] (mandatory)
    ch_long_reads_spades  // [val(meta), path(fastq)]                (mandatory)

    main:

    ch_versions = channel.empty()
    ch_assembled_contigs = channel.empty()

    if (!params.single_end && !params.skip_spadeshybrid) {
        ch_short_reads_spades_tmp = ch_short_reads_spades.map { meta, reads -> [meta.id, meta, reads] }

        ch_reads_spadeshybrid = ch_long_reads_spades
            .map { meta, reads -> [meta.id, meta, reads] }
            .combine(ch_short_reads_spades_tmp, by: 0)
            .map { _id, meta_long, long_reads, meta_short, short_reads ->
                if (meta_long.lr_platform == "OXFORD_NANOPORE" || meta_long.lr_platform == "OXFORD_NANOPORE_HQ") {
                    [meta_short, short_reads, [], long_reads]
                }
                else {
                    // For PacBio
                    [meta_short, short_reads, long_reads, []]
                }
            }

        METASPADESHYBRID(ch_reads_spadeshybrid, [], [])
        ch_versions = ch_versions.mix(METASPADESHYBRID.out.versions)

        ch_spadeshybrid_assemblies = METASPADESHYBRID.out.scaffolds.map { meta, assembly ->
            def meta_new = meta + [assembler: "SPAdesHybrid"]
            [meta_new, assembly]
        }
        ch_assembled_contigs = ch_assembled_contigs.mix(ch_spadeshybrid_assemblies)
    }

    emit:
    assembled_contigs = ch_assembled_contigs
    versions          = ch_versions
}
