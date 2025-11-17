/*
 * LONGREAD_PREPROCESSING: Preprocessing and QC for long reads
 */

include { NANOPLOT as NANOPLOT_RAW         } from '../../../modules/nf-core/nanoplot/main'
include { NANOPLOT as NANOPLOT_FILTERED    } from '../../../modules/nf-core/nanoplot/main'
include { NANOLYSE                         } from '../../../modules/nf-core/nanolyse/main'
include { PORECHOP_PORECHOP                } from '../../../modules/nf-core/porechop/porechop/main'
include { PORECHOP_ABI                     } from '../../../modules/nf-core/porechop/abi/main'
include { FILTLONG                         } from '../../../modules/nf-core/filtlong'
include { CHOPPER                          } from '../../../modules/nf-core/chopper'
include { NANOQ                            } from '../../../modules/nf-core/nanoq'
include { CAT_FASTQ as CAT_FASTQ_LONGREADS } from '../../../modules/nf-core/cat/fastq/main'

// include other subworkflows here
include { LONGREAD_HOSTREMOVAL             } from '../hostremoval_longread/main'

workflow LONGREAD_PREPROCESSING {
    take:
    ch_raw_long_reads // [val(meta), path(fastq)]                (mandatory)
    ch_short_reads    // [val(meta), path(fastq1), path(fastq2)]
    ch_lambda_db      // [val(meta), path(fasta)]
    ch_host_fasta     // [val(meta), path(fasta)]
    val_skip_qc       // val(boolean)

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    NANOPLOT_RAW(
        ch_raw_long_reads
    )
    ch_versions = ch_versions.mix(NANOPLOT_RAW.out.versions)

    if (!params.assembly_input) {
        if (!params.skip_adapter_trimming && !val_skip_qc) {
            if (params.longread_adaptertrimming_tool && params.longread_adaptertrimming_tool == 'porechop_abi') {
                PORECHOP_ABI(
                    ch_raw_long_reads,
                    [],
                )
                ch_versions = ch_versions.mix(PORECHOP_ABI.out.versions)
                ch_long_reads = PORECHOP_ABI.out.reads
                ch_multiqc_files = ch_multiqc_files.mix(PORECHOP_ABI.out.log)
            }
            else if (params.longread_adaptertrimming_tool == 'porechop') {
                PORECHOP_PORECHOP(
                    ch_raw_long_reads
                )
                ch_versions = ch_versions.mix(PORECHOP_PORECHOP.out.versions)
                ch_long_reads = PORECHOP_PORECHOP.out.reads
                ch_multiqc_files = ch_multiqc_files.mix(PORECHOP_PORECHOP.out.log)
            }
        }
        else {
            ch_long_reads = ch_raw_long_reads
        }

        if (!params.keep_lambda && params.longread_filtering_tool != 'chopper' && !val_skip_qc) {
            NANOLYSE(
                ch_long_reads,
                ch_lambda_db,
            )
            ch_versions = ch_versions.mix(NANOLYSE.out.versions)
            ch_long_reads = NANOLYSE.out.fastq
        }
        if (!params.skip_longread_filtering && !val_skip_qc) {
            if (params.longread_filtering_tool == 'filtlong') {
                // join long and short reads by sample name
                ch_short_reads_tmp = ch_short_reads.map { meta, sr -> [meta.id, sr] }

                ch_short_and_long_reads = ch_long_reads
                    .map { meta, lr -> [meta.id, meta, lr] }
                    .join(ch_short_reads_tmp, by: 0, remainder: true)
                    .filter { row -> row[1] != null }  // filter out samples with no long reads
                    .map { _id, meta_lr, lr, sr -> [meta_lr, sr ? sr : [], lr] }
                // should not occur for single-end, since SPAdes (hybrid) does not support single-end

                FILTLONG(
                    ch_short_and_long_reads
                )
                ch_versions = ch_versions.mix(FILTLONG.out.versions)
                ch_long_reads = FILTLONG.out.reads
                ch_multiqc_files = ch_multiqc_files.mix(FILTLONG.out.log)
            }
            else if (params.longread_filtering_tool == 'nanoq') {
                NANOQ(
                    ch_long_reads,
                    'fastq.gz',
                )
                ch_versions = ch_versions.mix(NANOQ.out.versions)
                ch_long_reads = NANOQ.out.reads
                ch_multiqc_files = ch_multiqc_files.mix(NANOQ.out.stats)
            }
            else if (params.longread_filtering_tool == 'chopper') {
                CHOPPER(
                    ch_long_reads,
                    ch_lambda_db,
                )
                ch_versions = ch_versions.mix(CHOPPER.out.versions)
                ch_long_reads = CHOPPER.out.fastq
            }
        }

        // host removal long reads
        if (params.host_fasta || params.host_genome) {
            LONGREAD_HOSTREMOVAL(
                ch_long_reads,
                ch_host_fasta,
            )
            ch_versions = ch_versions.mix(LONGREAD_HOSTREMOVAL.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(LONGREAD_HOSTREMOVAL.out.multiqc_files)
            ch_long_reads = LONGREAD_HOSTREMOVAL.out.reads
        }

        /**
         * Conditions for *not* running NANOPLOT_FILTERED:
         * - No host removal and skip_qc (params.skip_longread_qc)
         * - No host removal and *all* --keep_lambda, --skip_adapter_trimming, --skip_longread_filtering
         */
        if (!(val_skip_qc && !(params.host_fasta || params.host_genome))) {
            if (!(params.skip_adapter_trimming && params.skip_longread_filtering && params.keep_lambda && !(params.host_fasta || params.host_genome))) {
                NANOPLOT_FILTERED(ch_long_reads)
                ch_versions = ch_versions.mix(NANOPLOT_FILTERED.out.versions)
            }
        }

        // Run merging
        ch_long_reads_forcat = ch_long_reads
            .map { meta, reads ->
                def meta_new = meta - meta.subMap('run')
                meta_new.single_end = true
                [meta_new, reads]
            }
            .groupTuple()
            .branch { _meta, reads ->
                cat: reads.size() >= 2
                skip_cat: true
            }
        CAT_FASTQ_LONGREADS(ch_long_reads_forcat.cat.map { meta, reads -> [meta, reads.flatten()] })
        ch_versions = ch_versions.mix(CAT_FASTQ_LONGREADS.out.versions)

        ch_long_reads = CAT_FASTQ_LONGREADS.out.reads.mix(ch_long_reads_forcat.skip_cat.map { meta, reads -> [meta, reads[0]] })
    }
    else {
        ch_long_reads = ch_raw_long_reads
    }

    emit:
    long_reads    = ch_long_reads
    versions      = ch_versions
    multiqc_files = ch_multiqc_files
}
