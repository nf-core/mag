include { QUAST } from '../../modules/local/quast'

workflow ASSEMBLY_QC {
    take:
    assemblies

    main:
    ch_versions      = Channel.empty()
    ch_quast_multiqc = Channel.empty()

    if (!params.skip_quast){
        QUAST ( assemblies )
        ch_quast_multiqc = QUAST.out.qc
        ch_versions      = ch_versions.mix(QUAST.out.versions.first())
    }

    emit:
    quast_multiqc = ch_quast_multiqc
    versions      = ch_versions
}
