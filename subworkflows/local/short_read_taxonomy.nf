include { CENTRIFUGE_DB_PREPARATION } from '../../modules/local/centrifuge_db_preparation'
include { CENTRIFUGE                } from '../../modules/local/centrifuge'
include { KRAKEN2_DB_PREPARATION    } from '../../modules/local/kraken2_db_preparation'
include { KRAKEN2                   } from '../../modules/local/kraken2'
include { KRONA_DB                  } from '../../modules/local/krona_db'
include { KRONA                     } from '../../modules/local/krona'

workflow SHORT_READ_TAXONOMY {
    take:
    short_reads

    main:
    ch_versions = Channel.empty()

    if(params.centrifuge_db){
        ch_centrifuge_db_file = Channel
            .value(file( "${params.centrifuge_db}" ))
    } else {
        ch_centrifuge_db_file = Channel.empty()
    }

    if(params.kraken2_db){
        ch_kraken2_db_file = Channel
            .value(file( "${params.kraken2_db}" ))
    } else {
        ch_kraken2_db_file = Channel.empty()
    }

    CENTRIFUGE_DB_PREPARATION ( ch_centrifuge_db_file )
    CENTRIFUGE (
        short_reads,
        CENTRIFUGE_DB_PREPARATION.out.db
    )
    ch_versions = ch_versions.mix(CENTRIFUGE.out.versions.first())

    KRAKEN2_DB_PREPARATION (
        ch_kraken2_db_file
    )
    KRAKEN2 (
        short_reads,
        KRAKEN2_DB_PREPARATION.out.db
    )
    ch_versions = ch_versions.mix(KRAKEN2.out.versions.first())

    if (( params.centrifuge_db || params.kraken2_db ) && !params.skip_krona){
        KRONA_DB ()
        ch_tax_classifications = CENTRIFUGE.out.results_for_krona.mix(KRAKEN2.out.results_for_krona)
            . map { classifier, meta, report ->
                def meta_new = meta.clone()
                meta_new.classifier  = classifier
                [ meta_new, report ]
            }
        KRONA (
            ch_tax_classifications,
            KRONA_DB.out.db.collect()
        )
        ch_versions = ch_versions.mix(KRONA.out.versions.first())
    }

    emit:
    versions = ch_versions
}
