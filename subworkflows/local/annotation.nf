include { PROKKA   } from '../../modules/nf-core/prokka/main'
include { PRODIGAL } from '../../modules/nf-core/prodigal/main'

workflow ANNOTATION {
    take:
    bins_unbins
    assemblies

    main:
    ch_versions = Channel.empty()


    /*
    Prodigal: Predict proteins
    */
    if (!params.skip_prodigal){
        PRODIGAL (
            assemblies,
            'gff'
        )
        ch_versions = ch_versions.mix(PRODIGAL.out.versions.first())
    }

    /*
    * Prokka: Genome annotation
    */
    ch_bins_for_prokka = bins_unbins.transpose()
        .map { meta, bin ->
            def meta_new = meta.clone()
            meta_new.id  = bin.getBaseName()
            [ meta_new, bin ]
        }

    if (!params.skip_prokka){
        PROKKA (
            ch_bins_for_prokka,
            [],
            []
        )
        ch_versions = ch_versions.mix(PROKKA.out.versions.first())
    }

    emit:
    versions = ch_versions
}
