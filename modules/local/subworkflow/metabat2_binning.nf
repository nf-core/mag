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
    assemblies           // channel: [ val(assemlber), val(name), path(assembly) ]
    reads                // channel: [ val(meta), [ reads ] ]

    main:
    // build bowtie2 index for all assemblies
    BOWTIE2_INDEX_ASSEMBLY ( assemblies )

    // combine assemblies with sample reads for binning depending on specified mapping mode
    // add dummy to allow combining by index
    ch_reads_bowtie2 = reads.map{ meta, reads -> ["dummy", meta.id, meta.group, reads] }
    if (params.binning_map_mode == 'all'){
        ch_bowtie2_input = BOWTIE2_INDEX_ASSEMBLY.out.assembly_index
            .combine(ch_reads_bowtie2)            // combine assemblies with reads of all samples
            .map{ assembler, assembly_meta, assembly, index, dummy, reads_id, reads_group, reads -> [assembler, assembly_meta.id, assembly, index, reads_id, reads] }
    } else if (params.binning_map_mode == 'group'){
        ch_bowtie2_input = BOWTIE2_INDEX_ASSEMBLY.out.assembly_index
            .map { assembler, meta, assembly, index -> [assembler, meta.id, meta.group, assembly, index] }
            .combine(ch_reads_bowtie2, by: 2)     // combine assemblies with reads of samples from same group
            .map{ group, assembler, assembly_id, assembly, index, dummy, reads_id, reads -> [assembler, assembly_id, assembly, index, reads_id, reads] }
    } else {
        ch_bowtie2_input = BOWTIE2_INDEX_ASSEMBLY.out.assembly_index
            .map { assembler, meta, assembly, index -> [assembler, meta.id, meta.group, assembly, index] }
            .combine(ch_reads_bowtie2, by: 1)     // combine assemblies (not co-assembled) with reads from own sample
            .map{ name, assembler, assembly_group, assembly, index, dummy, reads_group, reads -> [assembler, name, assembly, index, name, reads] }
    }
    // TODO use already meta?
    BOWTIE2_ASSEMBLY ( ch_bowtie2_input )
    // group mappings for one assembly
    ch_grouped_mappings = BOWTIE2_ASSEMBLY.out.mappings
        .groupTuple(by:[0,1])
        .map { assembler, assembly_id, assembly, bams, bais -> [assembler, assembly_id, assembly[0], bams, bais] }     // multiple symlinks to the same assembly -> use first

    METABAT2 ( ch_grouped_mappings )

    emit:
    bins                     = METABAT2.out.bins
    bowtie2_assembly_multiqc = BOWTIE2_ASSEMBLY.out.log.map { assembler, assembly_id, reads_id, log -> if (assembly_id == reads_id) {return [ log ]} }
    bowtie2_version          = BOWTIE2_ASSEMBLY.out.version
    metabat2_version         = METABAT2.out.version
}