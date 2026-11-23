/*
 * Binning with MetaBAT2 and MaxBin2
 */
include { FASTA_BINNING_CONCOCT                                                                  } from '../../../subworkflows/nf-core/fasta_binning_concoct/main'
include { BINNING_METABINNER                                                                     } from '../../../subworkflows/local/binning_metabinner/main'

include { METABAT2_METABAT2                                                                      } from '../../../modules/nf-core/metabat2/metabat2/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS as METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS_SHORTREAD } from '../../../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS as METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS_LONGREAD  } from '../../../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'
include { MAXBIN2                                                                                } from '../../../modules/nf-core/maxbin2/main'
include { COMEBIN_RUNCOMEBIN                                                                     } from '../../../modules/nf-core/comebin/runcomebin/main'
include { SEMIBIN_SINGLEEASYBIN                                                                  } from '../../../modules/nf-core/semibin/singleeasybin/main'

include { SEQKIT_STATS                                                                           } from '../../../modules/nf-core/seqkit/stats/main'

include { CONVERT_DEPTHS                                                                         } from '../../../modules/local/mag_depths_convert/main'
include { SPLIT_FASTA                                                                            } from '../../../modules/local/split_fasta/main'

workflow BINNING {
    take:
    ch_assemblies // [val(meta), path(assembly), path(bams), path(bais)]
    val_bin_min_size // val(int)
    val_bin_max_size // val(int)

    main:

    ch_versions = channel.empty()
    ch_input_splitfasta = channel.empty()
    ch_publish_metabat2 = channel.empty()
    ch_publish_maxbin2 = channel.empty()
    ch_publish_concoct_stats = channel.empty()
    ch_publish_concoct_bins = channel.empty()
    ch_publish_comebin = channel.empty()
    ch_publish_metabinner = channel.empty()
    ch_publish_semibin = channel.empty()

    // generate coverage depths for each contig and branch by assembler type
    ch_summarizedepth_input = ch_assemblies
        .map { meta, _assembly, bams, bais ->
            [meta, bams, bais]
        }
        .branch { meta, _bams, _bais ->
            longread: meta.assembler in ['FLYE', 'METAMDBG']
            shortread: true
        }

    // Process each through appropriate module
    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS_LONGREAD(ch_summarizedepth_input.longread)
    ch_versions = ch_versions.mix(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS_LONGREAD.out.versions)

    METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS_SHORTREAD(ch_summarizedepth_input.shortread)
    ch_versions = ch_versions.mix(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS_SHORTREAD.out.versions)

    // Merge the outputs
    ch_combined_depths = METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS_LONGREAD.out.depth.mix(
        METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS_SHORTREAD.out.depth
    )
    ch_metabat_depths = ch_combined_depths.map { meta, depths ->
        def meta_new = meta + [binner: 'MetaBAT2']
        [meta_new, depths]
    }

    // combine depths back with assemblies
    ch_metabat2_input = ch_assemblies
        .map { meta, assembly, bams, bais ->
            def meta_new = meta + [binner: 'MetaBAT2']
            [meta_new, assembly, bams, bais]
        }
        .join(ch_metabat_depths, by: 0)
        .map { meta, assembly, _bams, _bais, depths ->
            [meta, assembly, depths]
        }

    // convert metabat2 depth files to maxbin2
    if (!params.skip_maxbin2) {
        CONVERT_DEPTHS(ch_metabat2_input)
        ch_versions = ch_versions.mix(CONVERT_DEPTHS.out.versions)

        ch_maxbin2_input = CONVERT_DEPTHS.out.output.map { meta, assembly, reads, depth ->
            def meta_new = meta + [binner: 'MaxBin2']
            [meta_new, assembly, reads, depth]
        }
    }

    // per-binner-output channel for length filtering via seqkit stats
    ch_bins_for_seqkit = channel.empty()

    // MetaBAT2
    if (!params.skip_metabat2) {
        METABAT2_METABAT2(ch_metabat2_input)

        ch_bins_for_seqkit = ch_bins_for_seqkit.mix(METABAT2_METABAT2.out.fasta)
        ch_input_splitfasta = ch_input_splitfasta.mix(METABAT2_METABAT2.out.unbinned)
        ch_publish_metabat2 = METABAT2_METABAT2.out.fasta
            .join(METABAT2_METABAT2.out.tooshort, remainder: true)
            .join(METABAT2_METABAT2.out.lowdepth, remainder: true)
    }

    // MaxBin2
    if (!params.skip_maxbin2) {
        MAXBIN2(ch_maxbin2_input)
        ch_versions = ch_versions.mix(MAXBIN2.out.versions)

        ch_bins_for_seqkit = ch_bins_for_seqkit.mix(MAXBIN2.out.binned_fastas)
        ch_input_splitfasta = ch_input_splitfasta.mix(MAXBIN2.out.unbinned_fasta)
        ch_publish_maxbin2 = MAXBIN2.out.binned_fastas
            .join(MAXBIN2.out.tooshort_fasta)
            .join(MAXBIN2.out.summary)
            .join(MAXBIN2.out.abundance, remainder: true)
    }

    // CONCOCT
    if (!params.skip_concoct) {

        // CONCOCT fails with fewer than two contigs; countFasta stops reading after the second one
        ch_concoct_input = ch_assemblies
            .filter { meta, assembly, _bams, _bais ->
                if (assembly.countFasta(limit: 2) < 2) {
                    log.warn("[nf-core/mag]: Fewer than two contigs in ${meta.assembler} assembly of '${meta.id}', skipping CONCOCT.")
                    return false
                }
                return true
            }
            .map { meta, bins, bams, bais ->
                def meta_new = meta + [binner: 'CONCOCT']
                [meta_new, bins, bams, bais]
            }
            .multiMap { meta, bins, bams, bais ->
                bins: [meta, bins]
                bams: [meta, bams, bais]
            }

        FASTA_BINNING_CONCOCT(ch_concoct_input.bins, ch_concoct_input.bams)

        ch_bins_for_seqkit = ch_bins_for_seqkit.mix(FASTA_BINNING_CONCOCT.out.bins)
        ch_publish_concoct_stats = FASTA_BINNING_CONCOCT.out.coverage_table.mix(
            FASTA_BINNING_CONCOCT.out.original_csv,
            FASTA_BINNING_CONCOCT.out.raw_clustering_csv,
            FASTA_BINNING_CONCOCT.out.pca_original,
            FASTA_BINNING_CONCOCT.out.pca_transformed,
            FASTA_BINNING_CONCOCT.out.args_txt,
            FASTA_BINNING_CONCOCT.out.log_txt,
            FASTA_BINNING_CONCOCT.out.cluster_table,
        )
        // one bin per element, as bins are published flattened out of their subdirectory
        ch_publish_concoct_bins = FASTA_BINNING_CONCOCT.out.bins.flatMap { meta, bins -> [bins].flatten().collect { bin -> [meta, bin] } }
    }

    // COMEBin
    if (!params.skip_comebin) {
        ch_comebin_input = ch_assemblies.map { meta, assembly, bams, _bais ->
            def meta_new = meta + [binner: 'COMEBin']
            [meta_new, assembly, bams]
        }

        COMEBIN_RUNCOMEBIN(ch_comebin_input)
        ch_versions = ch_versions.mix(COMEBIN_RUNCOMEBIN.out.versions)

        ch_bins_for_seqkit = ch_bins_for_seqkit.mix(COMEBIN_RUNCOMEBIN.out.bins)
        ch_publish_comebin = COMEBIN_RUNCOMEBIN.out.bins
            .join(COMEBIN_RUNCOMEBIN.out.tsv)
            .join(COMEBIN_RUNCOMEBIN.out.log)
            .join(COMEBIN_RUNCOMEBIN.out.embeddings)
            .join(COMEBIN_RUNCOMEBIN.out.covembeddings)
    }

    // MetaBinner
    if (!params.skip_metabinner) {
        BINNING_METABINNER(
            ch_metabat2_input.map { meta, assembly, depths ->
                def meta_new = meta + [binner: 'MetaBinner']
                [meta_new, assembly, depths]
            }
        )
        ch_versions = ch_versions.mix(BINNING_METABINNER.out.versions)

        ch_bins_for_seqkit = ch_bins_for_seqkit.mix(BINNING_METABINNER.out.bins)
        ch_input_splitfasta = ch_input_splitfasta.mix(BINNING_METABINNER.out.unbinned)
        ch_publish_metabinner = BINNING_METABINNER.out.publish
    }

    // SemiBin2
    if (!params.skip_semibin) {
        ch_semibin_input = ch_assemblies.map { meta, assembly, bams, _bais ->
            def meta_new = meta + [binner: 'SemiBin2'] + [sample_count: bams.size()]
            [meta_new, assembly, bams]
        }

        SEMIBIN_SINGLEEASYBIN(
            ch_semibin_input
        )
        ch_versions = ch_versions.mix(SEMIBIN_SINGLEEASYBIN.out.versions)

        // must remove the additional metadata because "workflow DEPTHS" channel combination is sensitive to any additional fields!
        ch_semibin_bins = SEMIBIN_SINGLEEASYBIN.out.output_fasta.map { meta, bins ->
            def meta_new = meta - meta.subMap('sample_count')
            [meta_new, bins]
        }

        ch_bins_for_seqkit = ch_bins_for_seqkit.mix(ch_semibin_bins)
        ch_publish_semibin = SEMIBIN_SINGLEEASYBIN.out.output_fasta.join(SEMIBIN_SINGLEEASYBIN.out.log, remainder: true)
    }

    // Performance note: grouping seqkit jobs by sample across binners locks
    // downstream bin QC until the slowest binner finishes.
    // extract max length of all entries in each bin, to allow filtering out of too small bins
    SEQKIT_STATS(ch_bins_for_seqkit)
    ch_versions = ch_versions.mix(SEQKIT_STATS.out.versions)

    //
    // Logic: Gather all the bin lengths, then check if the number of bins after length
    //        filtering is 0. Error if so, but only if we had bins to begin with.
    //
    SEQKIT_STATS.out.stats
        .splitCsv(sep: '\t', header: true, strip: true)
        .map { _meta, row -> row.sum_len.toInteger() }
        .collect()
        .ifEmpty([])
        .subscribe { lengths ->
            def n_bins = lengths.size()
            def n_filtered_bins = lengths
                .findAll { bin_size ->
                    bin_size >= val_bin_min_size && (val_bin_max_size ? bin_size <= val_bin_max_size : true)
                }
                .size()
            if (n_bins > 0 && n_filtered_bins == 0) {
                error(
                    "[nf-core/mag] ERROR: no bins passed the bin size filter specified between " + "--bin_min_size ${val_bin_min_size} and " + "--bin_max_size ${val_bin_max_size}. Please adjust parameters."
                )
            }
        }

    // filter out too-short/too-long bins by total length within each binner output;
    // joining per binner-sample keeps bin QC unblocked by the slowest binner
    ch_binning_results = ch_bins_for_seqkit
        .join(SEQKIT_STATS.out.stats)
        .map { meta, bins, stats ->
            def lengths = stats
                .splitCsv(sep: '\t', header: true, strip: true)
                .collectEntries { row -> [(row.file): row.sum_len.toInteger()] }
            def kept_bins = [bins]
                .flatten()
                .findAll { bin ->
                    def bin_total_length = lengths[bin.name]
                    bin_total_length >= val_bin_min_size && (val_bin_max_size ? bin_total_length <= val_bin_max_size : true)
                }
            [meta, kept_bins]
        }
        .filter { _meta, bins -> bins }

    // remove too-short contigs from unbinned contigs
    SPLIT_FASTA(ch_input_splitfasta)
    ch_versions = ch_versions.mix(SPLIT_FASTA.out.versions)

    emit:
    bins                  = ch_binning_results
    unbinned              = SPLIT_FASTA.out.unbinned
    metabat2depths        = ch_combined_depths
    versions              = ch_versions
    publish               = ch_publish_metabat2.map { meta, bins, tooshort, lowdepth -> [meta, [metabat2_bins: bins, metabat2_tooshort: tooshort, metabat2_lowdepth: lowdepth]] }
        .mix(
            ch_publish_maxbin2.map { meta, bins, tooshort, summary, abundance -> [meta, [maxbin2_bins: bins, maxbin2_tooshort: tooshort, maxbin2_summary: summary, maxbin2_abundance: abundance]] },
            ch_publish_concoct_stats.map { meta, f -> [meta, [concoct_stats: f]] },
            ch_publish_concoct_bins.map { meta, f -> [meta, [concoct_bin: f]] },
            ch_publish_comebin.map { meta, bins, tsv, log, embeddings, covembeddings -> [meta, [comebin_bins: bins, comebin_tsv: tsv, comebin_log: log, comebin_embeddings: embeddings, comebin_covembeddings: covembeddings]] },
            ch_publish_metabinner.map { meta, membership, log, bins, unbinned -> [meta, [metabinner_membership: membership, metabinner_log: log, metabinner_bins: bins, metabinner_unbinned: unbinned]] },
            ch_publish_semibin.map { meta, bins, log -> [meta, [semibin_bins: bins, semibin_log: log]] },
            SPLIT_FASTA.out.unbinned.join(SPLIT_FASTA.out.pooled, remainder: true).map { meta, unbinned, pooled -> [meta, [unbinned: unbinned, unbinned_pooled: pooled]] },
        ) // [val(meta), val(record)]
}
