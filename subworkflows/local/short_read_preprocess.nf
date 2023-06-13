include { FASTP                                               } from '../../modules/nf-core/fastp/main'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_PE                 } from '../../modules/nf-core/adapterremoval/main'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_SE                 } from '../../modules/nf-core/adapterremoval/main'
include { BOWTIE2_REMOVAL_BUILD as BOWTIE2_HOST_REMOVAL_BUILD } from '../../modules/local/bowtie2_removal_build'
include { BOWTIE2_REMOVAL_ALIGN as BOWTIE2_HOST_REMOVAL_ALIGN } from '../../modules/local/bowtie2_removal_align'
include { BOWTIE2_REMOVAL_BUILD as BOWTIE2_PHIX_REMOVAL_BUILD } from '../../modules/local/bowtie2_removal_build'
include { BOWTIE2_REMOVAL_ALIGN as BOWTIE2_PHIX_REMOVAL_ALIGN } from '../../modules/local/bowtie2_removal_align'
include { FASTQC as FASTQC_TRIMMED                            } from '../../modules/nf-core/fastqc/main'
include { SEQTK_MERGEPE                                       } from '../../modules/nf-core/seqtk/mergepe/main'
include { BBMAP_BBNORM                                        } from '../../modules/nf-core/bbmap/bbnorm/main'
include { CAT_FASTQ                                           } from '../../modules/nf-core/cat/fastq/main'

workflow SHORT_READ_PREPROCESS {
    take:
    raw_short_reads

    main:
    ch_versions    = Channel.empty()

    // Databases and references
    if ( params.host_genome ) {
        host_fasta = params.genomes[params.host_genome].fasta ?: false
        ch_host_fasta = Channel
            .value(file( "${host_fasta}" ))
        host_bowtie2index = params.genomes[params.host_genome].bowtie2 ?: false
        ch_host_bowtie2index = Channel
            .value(file( "${host_bowtie2index}/*" ))
    } else if ( params.host_fasta ) {
        ch_host_fasta = Channel
            .value(file( "${params.host_fasta}" ))
    } else {
        ch_host_fasta = Channel.empty()
    }

    if(!params.keep_phix) {
        ch_phix_db_file = Channel
            .value(file( "${params.phix_reference}" ))
    }

    // MultiQC outputs
    ch_multiqc_readprep             = Channel.empty()
    ch_bowtie2_removal_host_multiqc = Channel.empty()
    ch_fastqc_trimmed_multiqc       = Channel.empty()

    // Preprocessing
    if ( !params.skip_clipping ) {
        if ( params.clip_tool == 'fastp' ) {
            ch_clipmerge_out = FASTP (
                raw_short_reads,
                [],
                params.fastp_save_trimmed_fail,
                []
            )
            ch_short_reads_prepped = FASTP.out.reads
            ch_multiqc_readprep    = FASTP.out.json.collect{it[1]}.ifEmpty([])
            ch_versions            = ch_versions.mix(FASTP.out.versions.first())

        } else if ( params.clip_tool == 'adapterremoval' ) {
            // due to strange output file scheme in AR2, have to manually separate
            // SE/PE to allow correct pulling of reads after.
            ch_adapterremoval_in = raw_short_reads
                .branch {
                        single: it[0]['single_end']
                        paired: !it[0]['single_end']
                    }

            ADAPTERREMOVAL_PE ( ch_adapterremoval_in.paired, [] )
            ADAPTERREMOVAL_SE ( ch_adapterremoval_in.single, [] )

            ch_short_reads_prepped = ch_short_reads.mix(ADAPTERREMOVAL_SE.out.singles_truncated, ADAPTERREMOVAL_PE.out.paired_truncated)
            ch_multiqc_readprep = ch_multiqc_readprep.mix (
                ADAPTERREMOVAL_PE.out.settings.collect{it[1]}.ifEmpty([]),
                ADAPTERREMOVAL_SE.out.settings.collect{it[1]}.ifEmpty([])
                )
            ch_versions = ch_versions.mix(ADAPTERREMOVAL_PE.out.versions.first(), ADAPTERREMOVAL_SE.out.versions.first())
        }
    } else {
        ch_short_reads_prepped = raw_short_reads
    }

    if (params.host_fasta && !params.assembly_input){
        BOWTIE2_HOST_REMOVAL_BUILD (
            ch_host_fasta
        )
        ch_host_bowtie2index = BOWTIE2_HOST_REMOVAL_BUILD.out.index
    }

    if (params.host_fasta || params.host_genome){
        BOWTIE2_HOST_REMOVAL_ALIGN (
            ch_short_reads_prepped,
            ch_host_bowtie2index
        )
        ch_short_reads_hostremoved = BOWTIE2_HOST_REMOVAL_ALIGN.out.reads
        ch_bowtie2_removal_host_multiqc = BOWTIE2_HOST_REMOVAL_ALIGN.out.log
        ch_versions = ch_versions.mix(BOWTIE2_HOST_REMOVAL_ALIGN.out.versions.first())
    } else {
        ch_short_reads_hostremoved = ch_short_reads_prepped
    }

    if(!params.keep_phix) {
        BOWTIE2_PHIX_REMOVAL_BUILD (
            ch_phix_db_file
        )
        BOWTIE2_PHIX_REMOVAL_ALIGN (
            ch_short_reads_hostremoved,
            BOWTIE2_PHIX_REMOVAL_BUILD.out.index
        )
        ch_short_reads_phixremoved = BOWTIE2_PHIX_REMOVAL_ALIGN.out.reads
        ch_versions = ch_versions.mix(BOWTIE2_PHIX_REMOVAL_ALIGN.out.versions.first())
    } else {
        ch_short_reads_phixremoved = ch_short_reads_hostremoved
    }

    if (!(params.keep_phix && params.skip_clipping && !(params.host_genome || params.host_fasta))) {
        FASTQC_TRIMMED (
            ch_short_reads_phixremoved
        )
        ch_versions = ch_versions.mix(FASTQC_TRIMMED.out.versions)
    }

    // Run/Lane merging

    ch_short_reads_forcat = ch_short_reads_phixremoved
        .map {
            meta, reads ->
                def meta_new = meta.clone()
                meta_new.remove('run')
            [ meta_new, reads ]
        }
        .groupTuple()
        .branch {
            meta, reads ->
                cat:       ( meta.single_end && reads.size() == 1 ) || ( !meta.single_end && reads.size() >= 2 )
                skip_cat: true // Can skip merging if only single lanes
        }

    CAT_FASTQ ( ch_short_reads_forcat.cat.map { meta, reads -> [ meta, reads.flatten() ]} )

    ch_short_reads = Channel.empty()
    ch_short_reads = CAT_FASTQ.out.reads.mix( ch_short_reads_forcat.skip_cat ).map { meta, reads -> [ meta, reads.flatten() ]}
    ch_versions    = ch_versions.mix(CAT_FASTQ.out.versions.first())

    if ( params.bbnorm ) {
        if ( params.coassemble_group ) {
            // Interleave pairs, to be able to treat them as single ends when calling bbnorm. This prepares
            // for dropping the single_end parameter, but keeps assembly modules as they are, i.e. not
            // accepting a mix of single end and pairs.
            SEQTK_MERGEPE (
                ch_short_reads.filter { ! it[0].single_end }
            )
            ch_versions = ch_versions.mix(SEQTK_MERGEPE.out.versions.first())
            // Combine the interleaved pairs with any single end libraries. Set the meta.single_end to true (used by the bbnorm module).
                ch_bbnorm = SEQTK_MERGEPE.out.reads
                    .mix(ch_short_reads.filter { it[0].single_end })
                    .map { [ [ id: sprintf("group%s", it[0].group), group: it[0].group, single_end: true ], it[1] ] }
                    .groupTuple()
        } else {
            ch_bbnorm = ch_short_reads
        }
        BBMAP_BBNORM ( ch_bbnorm )
        ch_versions = ch_versions.mix(BBMAP_BBNORM.out.versions)
        ch_short_reads_assembly = BBMAP_BBNORM.out.fastq
    } else {
        ch_short_reads_assembly = ch_short_reads
    }

    emit:
    short_reads                  = ch_short_reads
    short_reads_assembly         = params.bbnorm ? ch_short_reads_bbnorm : ch_short_reads
    multiqc_readprep             = ch_multiqc_readprep
    multiqc_bowtie2_removal_host = ch_bowtie2_removal_host_multiqc
    multiqc_fastqc_trimmed       = ch_fastqc_trimmed_multiqc
    versions                     = ch_versions
}
