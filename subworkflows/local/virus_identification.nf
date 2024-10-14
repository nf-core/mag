/*
 * geNomad: Identification of mobile genetic elements
 */

include { GENOMAD_DOWNLOAD   } from '../../modules/nf-core/genomad/download/main'
include { GENOMAD_ENDTOEND   } from '../../modules/nf-core/genomad/endtoend/main'

workflow VIRUS_IDENTIFICATION {
    take:
    ch_assemblies   // [ [ meta] , fasta    ], input scaffolds (mandatory)
    ch_genomad_db   // [ db                 ], presupplied geNomad database (optional)

    main:
    ch_versions = Channel.empty()

    if ( params.genomad_db ) {
        ch_db_for_genomad = ch_genomad_db
    } else {
        ch_db_for_genomad = GENOMAD_DOWNLOAD( ).genomad_db
        ch_versions.mix( GENOMAD_DOWNLOAD.out.versions )
    }

    ch_identified_viruses = GENOMAD_ENDTOEND ( ch_assemblies, ch_db_for_genomad ).virus_fasta
    ch_versions.mix( GENOMAD_ENDTOEND.out.versions )

    emit:
    identified_viruses = ch_identified_viruses
    versions = ch_versions

}
