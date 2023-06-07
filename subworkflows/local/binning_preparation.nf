/*
 * Binning preparation with Bowtie2
 */

include { BOWTIE2_ASSEMBLY_BUILD    } from '../../modules/local/bowtie2_assembly_build'
include { BOWTIE2_ASSEMBLY_ALIGN    } from '../../modules/local/bowtie2_assembly_align'

workflow BINNING_PREPARATION {
    take:
    assemblies           // channel: [ val(meta), path(assembly) ]
    reads                // channel: [ val(meta), [ reads ] ]

    main:
    // build bowtie2 index for all assemblies
    BOWTIE2_ASSEMBLY_BUILD ( assemblies )
    reads.dump(tag: "binning_prep_reads")
    // combine assemblies with sample reads for binning depending on specified mapping mode
    if (params.binning_map_mode == 'all'){
        // combine assemblies with reads of all samples
        ch_bowtie2_input = BOWTIE2_ASSEMBLY_BUILD.out.assembly_index
            .dump(tag: "binning_map_mode_all_pre_combine")
            .combine(reads)
            .dump(tag: "binning_map_mode_all")
    } else if (params.binning_map_mode == 'group'){
        // combine assemblies with reads of samples from same group
        ch_reads_bowtie2 = reads.dump(tag: "binning_map_mode_group_pre_map").map{ meta, reads -> [ meta.group, meta, reads ] }.dump(tag: "binning_map_mode_group_post_map")
        // PRINTS POST_MAP HERE, DOESN'T PRINT COMBINE
        ch_bowtie2_input = BOWTIE2_ASSEMBLY_BUILD.out.assembly_index
            .dump(tag: "binning_map_mode_group_assembly_build_premap_combine")
            .map { meta, assembly, index -> [ meta.group, meta, assembly, index ] }
            .dump(tag: "binning_map_mode_group_assembly_build_precombine")
            .combine(ch_reads_bowtie2, by: 0)
            .dump(tag: "binning_map_mode_group_assembly_build_premap")
            .map { group, assembly_meta, assembly, index, reads_meta, reads -> [ assembly_meta, assembly, index, reads_meta, reads ] }
            .dump(tag: "binning_map_mode_group_assembly_build_postmap")

    } else {
        // combine assemblies (not co-assembled) with reads from own sample
        ch_reads_bowtie2 = reads.dump(tag: "ch_reads_bowtie2_premap").map{ meta, reads -> [ meta.id, meta, reads ] }
        ch_bowtie2_input = BOWTIE2_ASSEMBLY_BUILD.out.assembly_index
            .dump(tag: "index")
            .map { meta, assembly, index -> [ meta.id, meta, assembly, index ] }
            .combine(ch_reads_bowtie2, by: 0)
            .map { id, assembly_meta, assembly, index, reads_meta, reads -> [ assembly_meta, assembly, index, reads_meta, reads ] }
            .dump(tag: "binning_map_mode_none")

    }

    BOWTIE2_ASSEMBLY_ALIGN ( ch_bowtie2_input )
    // group mappings for one assembly
    ch_grouped_mappings = BOWTIE2_ASSEMBLY_ALIGN.out.mappings
        .groupTuple(by: 0)
        .map { meta, assembly, bams, bais -> [ meta, assembly.sort()[0], bams, bais ] }     // multiple symlinks to the same assembly -> use first of sorted list

    emit:
    bowtie2_assembly_multiqc = BOWTIE2_ASSEMBLY_ALIGN.out.log.map { assembly_meta, reads_meta, log -> if (assembly_meta.id == reads_meta.id) {return [ log ]} }
    bowtie2_version          = BOWTIE2_ASSEMBLY_ALIGN.out.versions
    grouped_mappings         = ch_grouped_mappings
}
