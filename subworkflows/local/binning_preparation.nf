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

    // combine assemblies with sample reads for binning depending on specified mapping mode
    if (params.binning_map_mode == 'all'){
        // combine assemblies with reads of all samples
        ch_bowtie2_input = BOWTIE2_ASSEMBLY_BUILD.out.assembly_index
            .combine(reads)
    } else if (params.binning_map_mode == 'group'){
        // combine assemblies with reads of samples from same group
        ch_reads_bowtie2 = reads.map{ meta, sample_reads -> [ meta.group, meta, sample_reads ] }
        ch_bowtie2_input = BOWTIE2_ASSEMBLY_BUILD.out.assembly_index
            .map { meta, assembly, index -> [ meta.group, meta, assembly, index ] }
            .combine(ch_reads_bowtie2, by: 0)
            .map { _group, assembly_meta, assembly, index, reads_meta, sample_reads ->
                [ assembly_meta, assembly, index, reads_meta, sample_reads ]
            }

    } else {
        // i.e. --binning_map_mode 'own'
        // combine assemblies (not co-assembled) with reads from own sample
        ch_reads_bowtie2 = reads.map{ meta, sample_reads -> [ meta.id, meta, sample_reads ] }
        ch_bowtie2_input = BOWTIE2_ASSEMBLY_BUILD.out.assembly_index
            .map { meta, assembly, index -> [ meta.id, meta, assembly, index ] }
            .combine(ch_reads_bowtie2, by: 0)
            .map { _id, assembly_meta, assembly, index, reads_meta, sample_reads ->
                [ assembly_meta, assembly, index, reads_meta, sample_reads ]
            }

    }

    BOWTIE2_ASSEMBLY_ALIGN ( ch_bowtie2_input )
    // group mappings for one assembly
    ch_grouped_mappings = BOWTIE2_ASSEMBLY_ALIGN.out.mappings
        .groupTuple(by: 0)
        .map { meta, assembly, bams, bais -> [ meta, assembly.sort()[0], bams, bais ] }     // multiple symlinks to the same assembly -> use first of sorted list

    emit:
    bowtie2_assembly_multiqc = BOWTIE2_ASSEMBLY_ALIGN.out.log.map { _assembly_meta, _reads_meta, log -> [ log ] }
    bowtie2_version          = BOWTIE2_ASSEMBLY_ALIGN.out.versions
    grouped_mappings         = ch_grouped_mappings
}
