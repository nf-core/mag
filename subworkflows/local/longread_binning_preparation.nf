include { MINIMAP2_INDEX as MINIMAP2_ASSEMBLY_INDEX                 } from '../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN as MINIMAP2_ASSEMBLY_ALIGN                 } from '../../modules/nf-core/minimap2/align/main'

workflow LONGREAD_BINNING_PREPARATION {
    take:
    assemblies // channel: [ val(meta), path(assembly) ]
    reads      // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions       = Channel.empty()

    MINIMAP2_ASSEMBLY_INDEX ( assemblies )
    MINIMAP2_ASSEMBLY_ALIGN ( reads, MINIMAP2_ASSEMBLY_INDEX.out.index, true, 'bai', false, false )
    ch_versions      = ch_versions.mix( MINIMAP2_ASSEMBLY_ALIGN.out.versions.first() )

    ch_minimap2_aligned = MINIMAP2_ASSEMBLY_ALIGN.out.bam
        .join(MINIMAP2_ASSEMBLY_ALIGN.out.index, by: 0)
    ch_minimap2_aligned.view()

    if (params.binning_map_mode == 'all'){
        ch_grouped_mappings = ch_minimap2_aligned
            .combine(assemblies)
            .map { meta_reads, bam, bai, meta_assembly, assembly -> [ meta_assembly, assembly, bam, bai ] }
            .groupTuple(by: 0)
            .map { meta, assembly, bams, bais -> [ meta, assembly.sort()[0], bams, bais ] }
    } else if (params.binning_map_mode == 'group') {
        ch_grouped_mappings = ch_minimap2_aligned
            .map{reads_meta, bams, bai -> [reads_meta.group, reads_meta, bams, bai]}
            .join(
                assemblies
                    .map{assemblies_meta, assembly -> [ assemblies_meta.group, assemblies_meta, assembly ] }
                , by: 0
            )
            .map { group, reads_meta, bam, bai, assemblies_meta, assembly -> [ assemblies_meta, assembly, bam, bai ] }
            .groupTuple(by: 0)
            .map { meta, assembly, bams, bais -> [ meta, assembly.sort()[0], bams, bais ] }
    } else {
        ch_grouped_mappings = ch_minimap2_aligned
            .map{reads_meta, bams, bai -> [reads_meta.id, reads_meta, bams, bai]}
            .join(
                assemblies
                    .map{assemblies_meta, assembly -> [ assemblies_meta.id, assemblies_meta, assembly ] }
                , by: 0
            )
            .map { id, reads_meta, bam, bai, assemblies_meta, assembly -> [ assemblies_meta, assembly, bam, bai ] }
            .groupTuple(by: 0)
            .map { meta, assembly, bams, bais -> [ meta, assembly.sort()[0], bams, bais ] }
    }

    emit:
    versions           = ch_versions
    grouped_mappings   = ch_grouped_mappings


    // group mappings for one assembly


}
