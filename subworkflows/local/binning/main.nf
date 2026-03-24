/*
 * Binning with MetaBAT2 and MaxBin2
 */
include { FASTA_BINNING_CONCOCT                                     } from '../../../subworkflows/nf-core/fasta_binning_concoct/main'
include { BINNING_METABINNER                                         } from '../../../subworkflows/local/binning_metabinner/main'

include { METABAT2_METABAT2                                          } from '../../../modules/nf-core/metabat2/metabat2/main'
include { COVERM_CONTIG as COVERM_CONTIG_SHORTREAD } from '../../../modules/nf-core/coverm/contig/main'
include { COVERM_CONTIG as COVERM_CONTIG_LONGREAD  } from '../../../modules/nf-core/coverm/contig/main'
include { MAXBIN2                                                    } from '../../../modules/nf-core/maxbin2/main'
include { COMEBIN_RUNCOMEBIN                                         } from '../../../modules/nf-core/comebin/runcomebin/main'
include { SEMIBIN_SINGLEEASYBIN                                      } from '../../../modules/nf-core/semibin/singleeasybin/main'

include { GUNZIP as GUNZIP_BINS                                      } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_UNBINS                                    } from '../../../modules/nf-core/gunzip/main'
include { SEQKIT_STATS                                               } from '../../../modules/nf-core/seqkit/stats/main'

include { CONVERT_DEPTHS                                             } from '../../../modules/local/mag_depths_convert/main'
include { ADJUST_MAXBIN2_EXT                                         } from '../../../modules/local/adjust_maxbin2_ext/main'
include { SPLIT_FASTA                                                } from '../../../modules/local/split_fasta/main'

workflow BINNING {
    take:
    ch_assemblies    // [val(meta), path(assembly), path(bams_or_reads), path(bais)]
    val_bin_min_size // val(int)
    val_bin_max_size // val(int)

    main:

    ch_versions = channel.empty()
    ch_input_splitfasta = channel.empty()

    // Validate: BAM-requiring binners cannot be used with CoverM native mappers
    // because no BAM files are produced in that path
    if (params.coverage_mapper != 'bowtie2') {
        if (!params.skip_concoct) {
            error("[nf-core/mag] CONCOCT requires BAM files. Please use --skip_concoct or set --coverage_mapper 'bowtie2'.")
        }
        if (!params.skip_comebin) {
            error("[nf-core/mag] COMEBin requires BAM files. Please use --skip_comebin or set --coverage_mapper 'bowtie2'.")
        }
        if (!params.skip_semibin) {
            error("[nf-core/mag] SemiBin2 requires BAM files. Please use --skip_semibin or set --coverage_mapper 'bowtie2'.")
        }
    }

    // Branch assemblies by assembler type for appropriate depth estimation
    ch_assemblies_branched = ch_assemblies
        .branch { meta, _assembly, _files, _indices ->
            longread:  meta.assembler in ['FLYE', 'METAMDBG']
            shortread: true
        }

    // Long reads always come as BAMs (from minimap2)
    ch_longread_depth = ch_assemblies_branched.longread
        .multiMap { meta, _assembly, bams, _bais ->
            reads:     [meta, bams]
            reference: [meta, []]
        }
    COVERM_CONTIG_LONGREAD(ch_longread_depth.reads, ch_longread_depth.reference, true, false)
    ch_longread_contig_depths = COVERM_CONTIG_LONGREAD.out.coverage

    // Short reads: use CoverM for depth calculation
    if (params.coverage_mapper == 'bowtie2') {
        ch_shortread_depth = ch_assemblies_branched.shortread
            .multiMap { meta, _assembly, bams, _bais ->
                reads:     [meta, bams]
                reference: [meta, []]
            }
        COVERM_CONTIG_SHORTREAD(ch_shortread_depth.reads, ch_shortread_depth.reference, true, false)
    }
    else {
        // CoverM native mapper: pass reads + assembly reference directly
        ch_shortread_depth = ch_assemblies_branched.shortread
            .multiMap { meta, assembly, reads, _bais ->
                reads:     [meta, reads]
                reference: [meta, assembly]
            }
        COVERM_CONTIG_SHORTREAD(ch_shortread_depth.reads, ch_shortread_depth.reference, false, false)
    }
    ch_shortread_contig_depths = COVERM_CONTIG_SHORTREAD.out.coverage

    // Merge depth outputs from short and long reads
    ch_combined_depths = ch_longread_contig_depths.mix(ch_shortread_contig_depths)
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

    // convert coverm depth files to maxbin2 abundance format
    if (!params.skip_maxbin2) {
        CONVERT_DEPTHS(ch_metabat2_input)
        ch_versions = ch_versions.mix(CONVERT_DEPTHS.out.versions)

        ch_maxbin2_input = CONVERT_DEPTHS.out.output.map { meta, assembly, reads, depth ->
            def meta_new = meta + [binner: 'MaxBin2']
            [meta_new, assembly, reads, depth]
        }
    }

    // main bins for decompressing for MAG_DEPTHS
    ch_bins_for_seqkit = channel.empty()

    // final gzipped bins
    ch_binning_results_gzipped_final = channel.empty()

    // MetaBAT2
    if (!params.skip_metabat2) {
        METABAT2_METABAT2(ch_metabat2_input)

        // before decompressing first have to separate and re-group due to limitation of GUNZIP module
        ch_bins_for_seqkit = ch_bins_for_seqkit.mix(METABAT2_METABAT2.out.fasta.transpose())
        ch_binning_results_gzipped_final = ch_binning_results_gzipped_final.mix(METABAT2_METABAT2.out.fasta)
        ch_input_splitfasta = ch_input_splitfasta.mix(METABAT2_METABAT2.out.unbinned)
    }

    // MaxBin2
    if (!params.skip_maxbin2) {
        MAXBIN2(ch_maxbin2_input)
        ch_versions = ch_versions.mix(MAXBIN2.out.versions)

        ADJUST_MAXBIN2_EXT(MAXBIN2.out.binned_fastas)
        ch_versions = ch_versions.mix(ADJUST_MAXBIN2_EXT.out.versions)

        ch_bins_for_seqkit = ch_bins_for_seqkit.mix(ADJUST_MAXBIN2_EXT.out.renamed_bins.transpose())
        ch_binning_results_gzipped_final = ch_binning_results_gzipped_final.mix(ADJUST_MAXBIN2_EXT.out.renamed_bins)
        ch_input_splitfasta = ch_input_splitfasta.mix(MAXBIN2.out.unbinned_fasta)
    }

    // CONCOCT (requires BAMs - only available with --coverm_mapper bowtie2)
    if (!params.skip_concoct) {

        ch_concoct_input = ch_assemblies
            .map { meta, bins, bams, bais ->
                def meta_new = meta + [binner: 'CONCOCT']
                [meta_new, bins, bams, bais]
            }
            .multiMap { meta, bins, bams, bais ->
                bins: [meta, bins]
                bams: [meta, bams, bais]
            }

        FASTA_BINNING_CONCOCT(ch_concoct_input.bins, ch_concoct_input.bams)
        ch_versions = ch_versions.mix(FASTA_BINNING_CONCOCT.out.versions)

        ch_bins_for_seqkit = ch_bins_for_seqkit.mix(FASTA_BINNING_CONCOCT.out.bins.transpose())
        ch_binning_results_gzipped_final = ch_binning_results_gzipped_final.mix(FASTA_BINNING_CONCOCT.out.bins)
    }

    // COMEBin (requires BAMs - only available with --coverm_mapper bowtie2)
    if (!params.skip_comebin) {
        ch_comebin_input = ch_assemblies.map { meta, assembly, bams, _bais ->
            def meta_new = meta + [binner: 'COMEBin']
            [meta_new, assembly, bams]
        }

        COMEBIN_RUNCOMEBIN(ch_comebin_input)
        ch_versions = ch_versions.mix(COMEBIN_RUNCOMEBIN.out.versions)

        ch_bins_for_seqkit = ch_bins_for_seqkit.mix(COMEBIN_RUNCOMEBIN.out.bins.transpose())
        ch_binning_results_gzipped_final = ch_binning_results_gzipped_final.mix(COMEBIN_RUNCOMEBIN.out.bins)
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

        ch_bins_for_seqkit = ch_bins_for_seqkit.mix(BINNING_METABINNER.out.bins.transpose())
        ch_binning_results_gzipped_final = ch_binning_results_gzipped_final.mix(BINNING_METABINNER.out.bins)
        ch_input_splitfasta = ch_input_splitfasta.mix(BINNING_METABINNER.out.unbinned)
    }

    // SemiBin2 (requires BAMs - only available with --coverm_mapper bowtie2)
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

        ch_bins_for_seqkit = ch_bins_for_seqkit.mix(ch_semibin_bins.transpose())
        ch_binning_results_gzipped_final = ch_binning_results_gzipped_final.mix(ch_semibin_bins)
    }

    // group bins into per-sample process and not flood clusters with thousands of seqkit jobs
    ch_bins_for_seqkitstats = ch_bins_for_seqkit
        .map { meta, bin ->
            [[id: meta.id], bin]
        }
        .groupTuple(by: 0)

    // extract max length of all entries in each bin, to allow filtering out of too small bins
    SEQKIT_STATS(ch_bins_for_seqkitstats)
    ch_versions = ch_versions.mix(SEQKIT_STATS.out.versions)

    //
    ch_seqkitstats_results = SEQKIT_STATS.out.stats
        .splitCsv(sep: '\t', header: true, strip: true)
        .map { _meta, row ->
            [[filename: row.file], [bin_total_length: row.sum_len.toInteger()]]
        }

    //
    // Logic: Gather all the bin lengths, then check if the number of bins after length
    //        filtering is 0. Error if so, but only if we had bins to begin with.
    //
    ch_seqkitstats_results
        .map { _meta, stats -> stats.bin_total_length }
        .collect()
        .ifEmpty([])
        .subscribe { stats ->
            def n_bins = stats.size()
            def n_filtered_bins = stats
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

    ch_final_bins_for_gunzip = ch_bins_for_seqkit
        .map { meta, bin ->
            [[filename: bin.name], meta, bin]
        }
        .join(ch_seqkitstats_results)
        .map { _key, meta, bin, stats ->
            [meta + stats, bin]
        }
        .filter { meta, _bin ->
            meta.bin_total_length >= val_bin_min_size && (val_bin_max_size ? meta.bin_total_length <= val_bin_max_size : true)
        }
        .map { meta, bin ->
            [meta.minus([bin_total_length: meta.bin_total_length]), bin]
        }

    // remove too-short contigs from unbinned contigs
    SPLIT_FASTA(ch_input_splitfasta)
    ch_versions = ch_versions.mix(SPLIT_FASTA.out.versions)

    // large unbinned contigs from SPLIT_FASTA for decompressing for MAG_DEPTHS,
    // first have to separate and re-group due to limitation of GUNZIP module
    ch_split_fasta_results_transposed = SPLIT_FASTA.out.unbinned.transpose()

    GUNZIP_BINS(ch_final_bins_for_gunzip)
    ch_versions = ch_versions.mix(GUNZIP_BINS.out.versions)
    ch_binning_results_gunzipped = GUNZIP_BINS.out.gunzip.groupTuple(by: 0)

    GUNZIP_UNBINS(ch_split_fasta_results_transposed)
    ch_versions = ch_versions.mix(GUNZIP_UNBINS.out.versions)
    ch_splitfasta_results_gunzipped = GUNZIP_UNBINS.out.gunzip.groupTuple(by: 0)

    emit:
    bins          = ch_binning_results_gunzipped
    bins_gz       = ch_binning_results_gzipped_final
    unbinned      = ch_splitfasta_results_gunzipped
    unbinned_gz   = SPLIT_FASTA.out.unbinned
    contig_depths = ch_combined_depths
    versions      = ch_versions
}
