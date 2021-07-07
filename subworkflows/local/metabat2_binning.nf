/*
 * Binning with MetaBAT2
 */

params.bowtie2_build_options      = [:]
params.bowtie2_align_options      = [:]
params.metabat2_options           = [:]
params.mag_depths_options         = [:]
params.mag_depths_plot_options    = [:]
params.mag_depths_summary_options = [:]

include { BOWTIE2_ASSEMBLY_BUILD    } from '../../modules/local/bowtie2_assembly_build'   addParams( options: params.bowtie2_build_options      )
include { BOWTIE2_ASSEMBLY_ALIGN    } from '../../modules/local/bowtie2_assembly_align'   addParams( options: params.bowtie2_align_options      )
include { METABAT2                  } from '../../modules/local/metabat2'                 addParams( options: params.metabat2_options           )
include { MAG_DEPTHS                } from '../../modules/local/mag_depths'               addParams( options: params.mag_depths_options         )
include { MAG_DEPTHS_PLOT           } from '../../modules/local/mag_depths_plot'          addParams( options: params.mag_depths_plot_options    )
include { MAG_DEPTHS_SUMMARY        } from '../../modules/local/mag_depths_summary'       addParams( options: params.mag_depths_summary_options )

workflow METABAT2_BINNING {
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
        ch_reads_bowtie2 = reads.map{ meta, reads -> [ meta.group, meta, reads ] }
        ch_bowtie2_input = BOWTIE2_ASSEMBLY_BUILD.out.assembly_index
            .map { meta, assembly, index -> [ meta.group, meta, assembly, index ] }
            .combine(ch_reads_bowtie2, by: 0)
            .map { group, assembly_meta, assembly, index, reads_meta, reads -> [ assembly_meta, assembly, index, reads_meta, reads ] }
    } else {
        // combine assemblies (not co-assembled) with reads from own sample
        ch_reads_bowtie2 = reads.map{ meta, reads -> [ meta.id, meta, reads ] }
        ch_bowtie2_input = BOWTIE2_ASSEMBLY_BUILD.out.assembly_index
            .map { meta, assembly, index -> [ meta.id, meta, assembly, index ] }
            .combine(ch_reads_bowtie2, by: 0)
            .map { id, assembly_meta, assembly, index, reads_meta, reads -> [ assembly_meta, assembly, index, reads_meta, reads ] }
    }

    BOWTIE2_ASSEMBLY_ALIGN ( ch_bowtie2_input )
    // group mappings for one assembly
    ch_grouped_mappings = BOWTIE2_ASSEMBLY_ALIGN.out.mappings
        .groupTuple(by: 0)
        .map { meta, assembly, bams, bais -> [ meta, assembly[0], bams, bais ] }     // multiple symlinks to the same assembly -> use first

    METABAT2 ( ch_grouped_mappings )

    // Compute bin depths for different samples (according to `binning_map_mode`)
    MAG_DEPTHS (
        METABAT2.out.bins,
        METABAT2.out.depths
    )
    // Plot bin depths heatmap for each assembly and mapped samples (according to `binning_map_mode`)
    ch_sample_groups = reads
        .collectFile(name:'sample_groups.tsv'){ meta, reads -> meta.id + '\t' + meta.group + '\n' }
    MAG_DEPTHS_PLOT (
        MAG_DEPTHS.out.depths,
        ch_sample_groups.collect()
    )

    MAG_DEPTHS_SUMMARY ( MAG_DEPTHS.out.depths.map{it[1]}.collect() )

    emit:
    bowtie2_assembly_multiqc = BOWTIE2_ASSEMBLY_ALIGN.out.log.map { assembly_meta, reads_meta, log -> if (assembly_meta.id == reads_meta.id) {return [ log ]} }
    bowtie2_version          = BOWTIE2_ASSEMBLY_ALIGN.out.version
    bins                     = METABAT2.out.bins
    depths_summary           = MAG_DEPTHS_SUMMARY.out.summary
    metabat2_version         = METABAT2.out.version
}
