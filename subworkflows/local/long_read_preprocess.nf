include { PORECHOP                      } from '../../modules/local/porechop'
include { NANOLYSE                      } from '../../modules/local/nanolyse'
include { FILTLONG                      } from '../../modules/local/filtlong'
include { NANOPLOT as NANOPLOT_RAW      } from '../../modules/local/nanoplot'
include { NANOPLOT as NANOPLOT_FILTERED } from '../../modules/local/nanoplot'

workflow LONG_READ_PREPROCESS {
    take:
    raw_long_reads
    short_reads

    main:
    ch_versions = Channel.empty()

    // Databases and references
    if (!params.keep_lambda) {
        ch_nanolyse_db = Channel
            .value(file( "${params.lambda_reference}" ))
    }

    // Preprocessing
    NANOPLOT_RAW (
        raw_long_reads
    )
    ch_versions = ch_versions.mix(NANOPLOT_RAW.out.versions.first())

    ch_long_reads = raw_long_reads

    if (!params.skip_adapter_trimming) {
        PORECHOP (
            ch_long_reads
        )
        ch_long_reads = PORECHOP.out.reads
        ch_versions = ch_versions.mix(PORECHOP.out.versions.first())
    }

    if (!params.keep_lambda) {
        NANOLYSE (
            ch_long_reads,
            ch_nanolyse_db
        )
        ch_long_reads = NANOLYSE.out.reads
        ch_versions = ch_versions.mix(NANOLYSE.out.versions.first())
    }

    // join long and short reads by sample name
    ch_short_reads_tmp = short_reads
        .map { meta, sr -> [ meta.id, meta, sr ] }

    ch_short_and_long_reads = ch_long_reads
        .map { meta, lr -> [ meta.id, meta, lr ] }
        .join(ch_short_reads_tmp, by: 0)
        .map { id, meta_lr, lr, meta_sr, sr -> [ meta_lr, lr, sr[0], sr[1] ] }  // should not occur for single-end, since SPAdes (hybrid) does not support single-end

    FILTLONG (
        ch_short_and_long_reads
    )
    ch_long_reads = FILTLONG.out.reads
    ch_versions = ch_versions.mix(FILTLONG.out.versions.first())

    NANOPLOT_FILTERED (
        ch_long_reads
    )

    emit:
    long_reads = ch_long_reads
    versions   = ch_versions
}
