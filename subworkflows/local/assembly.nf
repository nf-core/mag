include { POOL_SINGLE_READS as POOL_SHORT_SINGLE_READS        } from '../../modules/local/pool_single_reads'
include { POOL_PAIRED_READS                                   } from '../../modules/local/pool_paired_reads'
include { POOL_SINGLE_READS as POOL_LONG_READS                } from '../../modules/local/pool_single_reads'
include { MEGAHIT                                             } from '../../modules/local/megahit'
include { SPADES                                              } from '../../modules/local/spades'
include { SPADESHYBRID                                        } from '../../modules/local/spadeshybrid'

workflow ASSEMBLY {
    take:
    short_reads_assembly
    long_reads

    main:
    ch_versions = Channel.empty()

    // Co-assembly: prepare grouping for MEGAHIT and for pooling for SPAdes
    if (params.coassemble_group) {
        // short reads
        // group and set group as new id
        ch_short_reads_grouped = short_reads_assembly
            .map { meta, reads -> [ meta.group, meta, reads ] }
            .groupTuple(by: 0)
            .map { group, metas, reads ->
                def assemble_as_single = params.single_end || ( params.bbnorm && params.coassemble_group )
                def meta         = [:]
                meta.id          = "group-$group"
                meta.group       = group
                meta.single_end  = assemble_as_single
                if ( assemble_as_single ) [ meta, reads.collect { it }, [] ]
                else [ meta, reads.collect { it[0] }, reads.collect { it[1] } ]
            }
        // long reads
        // group and set group as new id
        ch_long_reads_grouped = long_reads
            .map { meta, reads -> [ meta.group, meta, reads ] }
            .groupTuple(by: 0)
            .map { group, metas, reads ->
                def meta = [:]
                meta.id          = "group-$group"
                meta.group       = group
                [ meta, reads.collect { it } ]
            }
    } else {
        ch_short_reads_grouped = short_reads_assembly
            .filter { it[0].single_end }
            .map { meta, reads -> [ meta, [ reads ], [] ] }
            .mix (
                short_reads_assembly
                    .filter { ! it[0].single_end }
                    .map { meta, reads -> [ meta, [ reads[0] ], [ reads[1] ] ] }
            )
        ch_long_reads_grouped = long_reads
    }

    ch_assemblies = Channel.empty()
    if (!params.skip_megahit){
        MEGAHIT ( ch_short_reads_grouped )
        ch_megahit_assemblies = MEGAHIT.out.assembly
            .map { meta, assembly ->
                def meta_new = meta.clone()
                meta_new.assembler  = "MEGAHIT"
                [ meta_new, assembly ]
            }
        ch_assemblies = ch_assemblies.mix(ch_megahit_assemblies)
        ch_versions = ch_versions.mix(MEGAHIT.out.versions.first())
    }

    // Co-assembly: pool reads for SPAdes
    if ( ! params.skip_spades || ! params.skip_spadeshybrid ){
        if ( params.coassemble_group ) {
            if ( params.bbnorm ) {
                ch_short_reads_spades = ch_short_reads_grouped.map { [ it[0], it[1] ] }
            } else {
                POOL_SHORT_SINGLE_READS (
                    ch_short_reads_grouped
                        .filter { it[0].single_end }
                )
                POOL_PAIRED_READS (
                    ch_short_reads_grouped
                        .filter { ! it[0].single_end }
                )
                ch_short_reads_spades = POOL_SHORT_SINGLE_READS.out.reads
                    .mix(POOL_PAIRED_READS.out.reads)
            }
        } else {
            ch_short_reads_spades = short_reads_assembly
        }
        // long reads
        if (!params.single_end && !params.skip_spadeshybrid){
            POOL_LONG_READS ( ch_long_reads_grouped )
            ch_long_reads_spades = POOL_LONG_READS.out.reads
        } else {
            ch_long_reads_spades = Channel.empty()
        }
    } else {
        ch_short_reads_spades = Channel.empty()
        ch_long_reads_spades  = Channel.empty()
    }

    if (!params.single_end && !params.skip_spades){
        SPADES ( ch_short_reads_spades )
        ch_spades_assemblies = SPADES.out.assembly
            .map { meta, assembly ->
                def meta_new = meta.clone()
                meta_new.assembler  = "SPAdes"
                [ meta_new, assembly ]
            }
        ch_assemblies = ch_assemblies.mix(ch_spades_assemblies)
        ch_versions = ch_versions.mix(SPADES.out.versions.first())
    }

    if (!params.single_end && !params.skip_spadeshybrid){
        ch_short_reads_spades_tmp = ch_short_reads_spades
            .map { meta, reads -> [ meta.id, meta, reads ] }
        ch_reads_spadeshybrid = ch_long_reads_spades
            .map { meta, reads -> [ meta.id, meta, reads ] }
            .combine(ch_short_reads_spades_tmp, by: 0)
            .map { id, meta_long, long_reads, meta_short, short_reads -> [ meta_short, long_reads, short_reads ] }
        SPADESHYBRID ( ch_reads_spadeshybrid )
        ch_spadeshybrid_assemblies = SPADESHYBRID.out.assembly
            .map { meta, assembly ->
                def meta_new = meta.clone()
                meta_new.assembler  = "SPAdesHybrid"
                [ meta_new, assembly ]
            }
        ch_assemblies = ch_assemblies.mix(ch_spadeshybrid_assemblies)
        ch_versions = ch_versions.mix(SPADESHYBRID.out.versions.first())
    }

    emit:
    assemblies = ch_assemblies
    versions   = ch_versions
}
