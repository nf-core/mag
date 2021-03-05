/*
 * Binning with MetaBAT2
 */

params.bowtie2_index_options = [:]
params.bowtie2_align_options = [:]
params.metabat2_options      = [:]

include { BOWTIE2_INDEX_ASSEMBLY    } from '../process/bowtie2_index_assembly'   addParams( options: params.bowtie2_index_options )
include { BOWTIE2_ASSEMBLY          } from '../process/bowtie2_assembly'         addParams( options: params.bowtie2_align_options )
include { METABAT2                  } from '../process/metabat2'                 addParams( options: params.metabat2_options      )

workflow METABAT2_BINNING {
    take:
    assemblies           // channel: [ val(meta), path(assembly) ]
    reads                // channel: [ val(meta), [ reads ] ]

    main:
    // build bowtie2 index for all assemblies
    BOWTIE2_INDEX_ASSEMBLY ( assemblies )

    // combine assemblies with sample reads for binning depending on specified mapping mode
    if (params.binning_map_mode == 'all'){
        // combine assemblies with reads of all samples
        ch_bowtie2_input = BOWTIE2_INDEX_ASSEMBLY.out.assembly_index
            .combine(reads)
    } else if (params.binning_map_mode == 'group'){
        // combine assemblies with reads of samples from same group
        ch_reads_bowtie2 = reads.map{ meta, reads -> [ meta.group, meta, reads ] }
        ch_bowtie2_input = BOWTIE2_INDEX_ASSEMBLY.out.assembly_index
            .map { meta, assembly, index -> [ meta.group, meta, assembly, index ] }
            .combine(ch_reads_bowtie2, by: 0)
            .map { group, assembly_meta, assembly, index, reads_meta, reads -> [ assembly_meta, assembly, index, reads_meta, reads ] }
    } else {
        // combine assemblies (not co-assembled) with reads from own sample
        ch_reads_bowtie2 = reads.map{ meta, reads -> [ meta.id, meta, reads ] }
        ch_bowtie2_input = BOWTIE2_INDEX_ASSEMBLY.out.assembly_index
            .map { meta, assembly, index -> [ meta.id, meta, assembly, index ] }
            .combine(ch_reads_bowtie2, by: 0)
            .map { id, assembly_meta, assembly, index, reads_meta, reads -> [ assembly_meta, assembly, index, reads_meta, reads ] }
    }

    BOWTIE2_ASSEMBLY ( ch_bowtie2_input )
    // group mappings for one assembly
    ch_grouped_mappings = BOWTIE2_ASSEMBLY.out.mappings
        .groupTuple(by: 0)
        .map { meta, assembly, bams, bais -> [ meta, assembly[0], bams, bais ] }     // multiple symlinks to the same assembly -> use first

    METABAT2 ( ch_grouped_mappings )

    emit:
    bowtie2_assembly_multiqc = BOWTIE2_ASSEMBLY.out.log.map { assembly_meta, reads_meta, log -> if (assembly_meta.id == reads_meta.id) {return [ log ]} }
    bowtie2_version          = BOWTIE2_ASSEMBLY.out.version
    bins                     = METABAT2.out.bins
    metabat2_version         = METABAT2.out.version
}
