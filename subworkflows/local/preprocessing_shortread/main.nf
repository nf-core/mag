/*
 * SHORTREAD_PREPROCESSING: Preprocessing and QC for short reads
 */

include { FASTQC as FASTQC_RAW                                } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED                            } from '../../../modules/nf-core/fastqc/main'
include { FASTP                                               } from '../../../modules/nf-core/fastp/main'
include { TRIMMOMATIC                                         } from '../../../modules/nf-core/trimmomatic/main'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_PE                 } from '../../../modules/nf-core/adapterremoval/main'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_SE                 } from '../../../modules/nf-core/adapterremoval/main'
include { CAT_FASTQ                                           } from '../../../modules/nf-core/cat/fastq/main'
include { SEQTK_MERGEPE                                       } from '../../../modules/nf-core/seqtk/mergepe/main'
include { BBMAP_BBNORM                                        } from '../../../modules/nf-core/bbmap/bbnorm/main'

include { BOWTIE2_REMOVAL_BUILD as BOWTIE2_HOST_REMOVAL_BUILD } from '../../../modules/local/bowtie2_removal_build/main'
include { BOWTIE2_REMOVAL_ALIGN as BOWTIE2_HOST_REMOVAL_ALIGN } from '../../../modules/local/bowtie2_removal_align/main'
include { BOWTIE2_REMOVAL_BUILD as BOWTIE2_PHIX_REMOVAL_BUILD } from '../../../modules/local/bowtie2_removal_build/main'
include { BOWTIE2_REMOVAL_ALIGN as BOWTIE2_PHIX_REMOVAL_ALIGN } from '../../../modules/local/bowtie2_removal_align/main'

workflow SHORTREAD_PREPROCESSING {
    take:
    ch_raw_short_reads   // [val(meta), [path(fastq1), path(fastq2)]] (mandatory)
    ch_host_fasta        // [val(meta), path(fasta)]                  (optional)
    ch_host_genome_index // path(fasta)                               (optional)
    ch_phix_db_file      // [val(meta), path(fasta)]                  (optional)
    val_skip_qc          // val(boolean)
    val_skip_fastqc      // val(boolean)

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_publish_fastqc = channel.empty()
    ch_publish_clipping = channel.empty()
    ch_publish_clipping_reads = channel.empty()
    ch_publish_hostremoval = channel.empty()
    ch_publish_hostremoval_reads = channel.empty()
    ch_publish_phixremoval = channel.empty()
    ch_publish_phixremoval_reads = channel.empty()
    ch_publish_bbnorm = channel.empty()

    if (!val_skip_fastqc) {
        FASTQC_RAW(ch_raw_short_reads)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip)
        ch_publish_fastqc = ch_publish_fastqc.mix(FASTQC_RAW.out.html)
    }

    if (!params.assembly_input) {
        if (!params.skip_clipping && !val_skip_qc) {
            if (params.clip_tool == 'fastp') {
                FASTP(
                    ch_raw_short_reads.map { meta, reads -> [meta, reads, []] },
                    false,
                    params.fastp_save_trimmed_fail,
                    false,
                )
                ch_versions = ch_versions.mix(FASTP.out.versions)
                ch_short_reads_prepped = FASTP.out.reads
                ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json)
                ch_publish_clipping = ch_publish_clipping.mix(FASTP.out.html, FASTP.out.json)
                ch_publish_clipping_reads = ch_publish_clipping_reads.mix(FASTP.out.reads, FASTP.out.reads_fail)
            }
            else if (params.clip_tool == 'adapterremoval') {

                // due to strange output file scheme in AR2, have to manually separate
                // SE/PE to allow correct pulling of reads after.
                ch_adapterremoval_in = ch_raw_short_reads.branch { meta, _reads ->
                    single: meta.single_end
                    paired: !meta.single_end
                }

                ADAPTERREMOVAL_PE(ch_adapterremoval_in.paired, [])
                ch_versions = ch_versions.mix(ADAPTERREMOVAL_PE.out.versions)
                ADAPTERREMOVAL_SE(ch_adapterremoval_in.single, [])
                ch_versions = ch_versions.mix(ADAPTERREMOVAL_SE.out.versions)

                ch_short_reads_prepped = channel.empty()
                ch_short_reads_prepped = ch_short_reads_prepped.mix(
                    ADAPTERREMOVAL_SE.out.singles_truncated,
                    ADAPTERREMOVAL_PE.out.paired_truncated,
                )

                ch_multiqc_files = ch_multiqc_files.mix(ADAPTERREMOVAL_PE.out.settings)
                ch_multiqc_files = ch_multiqc_files.mix(ADAPTERREMOVAL_SE.out.settings)
                ch_publish_clipping = ch_publish_clipping.mix(ADAPTERREMOVAL_PE.out.settings, ADAPTERREMOVAL_SE.out.settings)
                ch_publish_clipping_reads = ch_publish_clipping_reads.mix(
                    ADAPTERREMOVAL_PE.out.paired_truncated,
                    ADAPTERREMOVAL_PE.out.singles_truncated,
                    ADAPTERREMOVAL_PE.out.collapsed_truncated,
                    ADAPTERREMOVAL_PE.out.discarded,
                )
            }
            else if (params.clip_tool == 'trimmomatic') {

                TRIMMOMATIC(ch_raw_short_reads)
                ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions)

                ch_short_reads_prepped = channel.empty()
                ch_short_reads_prepped = TRIMMOMATIC.out.trimmed_reads

                ch_multiqc_files = ch_multiqc_files.mix(TRIMMOMATIC.out.out_log)
                ch_publish_clipping = ch_publish_clipping.mix(TRIMMOMATIC.out.trim_log)
                ch_publish_clipping_reads = ch_publish_clipping_reads.mix(TRIMMOMATIC.out.trimmed_reads, TRIMMOMATIC.out.unpaired_reads)
            }
        }
        else {
            ch_short_reads_prepped = ch_raw_short_reads
        }

        if (params.host_fasta || params.host_genome) {
            if (params.host_fasta_bowtie2index) {
                ch_host_bowtie2index = file(params.host_fasta_bowtie2index, checkIfExists: true)
            }
            else {
                ch_host_fasta_for_build = ch_host_fasta
                    .combine(ch_short_reads_prepped)
                    .map { host_fasta, _meta, _reads ->
                        host_fasta
                    }
                    .first()
                // makes sure to only use the host fasta if the short read channel is not empty
                BOWTIE2_HOST_REMOVAL_BUILD(
                    ch_host_fasta_for_build
                )
                ch_versions = ch_versions.mix(BOWTIE2_HOST_REMOVAL_BUILD.out.versions)

                ch_host_bowtie2index = BOWTIE2_HOST_REMOVAL_BUILD.out.index
            }
        }
        else if (params.host_genome) {
            ch_host_bowtie2index = ch_host_genome_index
        }

        if (params.host_fasta || params.host_genome) {
            BOWTIE2_HOST_REMOVAL_ALIGN(
                ch_short_reads_prepped,
                ch_host_bowtie2index,
            )
            ch_versions = ch_versions.mix(BOWTIE2_HOST_REMOVAL_ALIGN.out.versions)

            ch_short_reads_hostremoved = BOWTIE2_HOST_REMOVAL_ALIGN.out.reads
            ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_HOST_REMOVAL_ALIGN.out.log)
            ch_publish_hostremoval = ch_publish_hostremoval.mix(BOWTIE2_HOST_REMOVAL_ALIGN.out.log, BOWTIE2_HOST_REMOVAL_ALIGN.out.read_ids)
            ch_publish_hostremoval_reads = BOWTIE2_HOST_REMOVAL_ALIGN.out.reads
        }
        else {
            ch_short_reads_hostremoved = ch_short_reads_prepped
        }

        if (!params.keep_phix && !val_skip_qc) {
            ch_phix_fasta_for_build = ch_phix_db_file
                .combine(ch_short_reads_prepped)
                .map { host_fasta, _meta, _reads ->
                    host_fasta
                }
                .first()
            BOWTIE2_PHIX_REMOVAL_BUILD(
                ch_phix_fasta_for_build
            )
            ch_versions = ch_versions.mix(BOWTIE2_PHIX_REMOVAL_BUILD.out.versions)

            BOWTIE2_PHIX_REMOVAL_ALIGN(
                ch_short_reads_hostremoved,
                BOWTIE2_PHIX_REMOVAL_BUILD.out.index,
            )
            ch_versions = ch_versions.mix(BOWTIE2_PHIX_REMOVAL_ALIGN.out.versions)

            ch_short_reads_for_merge = BOWTIE2_PHIX_REMOVAL_ALIGN.out.reads
            ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_PHIX_REMOVAL_ALIGN.out.log)
            ch_publish_phixremoval = ch_publish_phixremoval.mix(BOWTIE2_PHIX_REMOVAL_ALIGN.out.log)
            ch_publish_phixremoval_reads = BOWTIE2_PHIX_REMOVAL_ALIGN.out.reads
        }
        else {
            ch_short_reads_for_merge = ch_short_reads_hostremoved
        }

        /**
         * Conditions for *not* running FASTQC_TRIMMED:
         * - Skip fastqc (params.skip_fastqc)
         * - No host removal and skip_qc (params.skip_shortread_qc)
         * - No host removal and *both* --keep_phix --skip_clipping
         */
        if (!val_skip_fastqc) {
            if (!(val_skip_qc && !(params.host_fasta || params.host_genome))) {
                if (!(params.keep_phix && params.skip_clipping && !(params.host_genome || params.host_fasta))) {
                    FASTQC_TRIMMED(
                        ch_short_reads_for_merge
                    )
                    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip)
                    ch_publish_fastqc = ch_publish_fastqc.mix(FASTQC_TRIMMED.out.html)
                }
            }
        }
    }
    else {
        ch_short_reads_for_merge = ch_raw_short_reads
    }

    // Run/Lane merging

    ch_short_reads_forcat = ch_short_reads_for_merge
        .map { meta, reads ->
            def meta_new = meta - meta.subMap('run')
            [meta_new, reads]
        }
        .groupTuple()
        .branch { _meta, reads ->
            cat: reads.size() >= 2
            skip_cat: true
        }

    CAT_FASTQ(ch_short_reads_forcat.cat.map { meta, reads -> [meta, reads.flatten()] })
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)


    // Ensure we don't have nests of nests so that structure is in form expected for assembly
    ch_short_reads_catskipped = ch_short_reads_forcat.skip_cat.map { meta, reads ->
        def new_reads = meta.single_end ? reads[0] : reads.flatten()
        [meta, new_reads]
    }

    // Combine single run and multi-run-merged data
    ch_short_reads = CAT_FASTQ.out.reads.mix(ch_short_reads_catskipped)

    if (params.bbnorm && !params.assembly_input) {
        if (params.coassemble_group) {
            // Interleave pairs, to be able to treat them as single ends when calling bbnorm. This prepares
            // for dropping the single_end parameter, but keeps assembly modules as they are, i.e. not
            // accepting a mix of single end and pairs.
            SEQTK_MERGEPE(
                ch_short_reads.filter { meta, _reads -> !meta.single_end }
            )
            ch_versions = ch_versions.mix(SEQTK_MERGEPE.out.versions)
            // Combine the interleaved pairs with any single end libraries. Set the meta.single_end to true (used by the bbnorm module).
            ch_bbnorm = SEQTK_MERGEPE.out.reads
                .mix(ch_short_reads.filter { meta, _reads -> meta.single_end })
                .map { meta, reads ->
                    [[id: "group${meta.group}", group: meta.group, single_end: true], reads]
                }
                .groupTuple()
        }
        else {
            ch_bbnorm = ch_short_reads
        }
        BBMAP_BBNORM(ch_bbnorm)
        ch_versions = ch_versions.mix(BBMAP_BBNORM.out.versions)
        ch_short_reads_assembly = BBMAP_BBNORM.out.fastq
        ch_publish_bbnorm = BBMAP_BBNORM.out.fastq.join(BBMAP_BBNORM.out.log)
    }
    else {
        ch_short_reads_assembly = ch_short_reads
    }

    emit:
    short_reads               = ch_short_reads
    short_reads_assembly      = ch_short_reads_assembly
    versions                  = ch_versions
    multiqc_files             = ch_multiqc_files
    publish                   = ch_publish_fastqc.map { meta, f -> [meta, [fastqc: f]] }
        .mix(
            ch_publish_clipping.map { meta, f -> [meta, [clipping: f]] },
            ch_publish_clipping_reads.map { meta, f -> [meta, [clipping_reads: f]] },
            ch_publish_hostremoval.map { meta, f -> [meta, [hostremoval: f]] },
            ch_publish_hostremoval_reads.map { meta, f -> [meta, [hostremoval_reads: f]] },
            ch_publish_phixremoval.map { meta, f -> [meta, [phixremoval: f]] },
            ch_publish_phixremoval_reads.map { meta, f -> [meta, [phixremoval_reads: f]] },
            ch_publish_bbnorm.map { meta, fastq, log -> [meta, [bbnorm_fastq: fastq, bbnorm_log: log]] },
        ) // [val(meta), val(record)]
}
