include { SHORTREAD_BINNING_PREPARATION } from '../binning_preparation_shortread/main'
include { LONGREAD_BINNING_PREPARATION  } from '../binning_preparation_longread/main'

workflow BINNING_PREPARATION {
    take:
    ch_shortread_assemblies // [val(meta), path(assembly)]
    ch_shortreads           // [val(meta), path(reads)]
    ch_longread_assemblies  // [val(meta), path(assembly)]
    ch_longreads            // [val(meta), path(reads)]

    main:
    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    SHORTREAD_BINNING_PREPARATION(ch_shortread_assemblies, ch_shortreads)
    ch_versions = ch_versions.mix(SHORTREAD_BINNING_PREPARATION.out.versions)

    LONGREAD_BINNING_PREPARATION(ch_longread_assemblies, ch_longreads)
    ch_versions = ch_versions.mix(LONGREAD_BINNING_PREPARATION.out.versions)

    ch_grouped_mappings = SHORTREAD_BINNING_PREPARATION.out.grouped_mappings.mix(
        LONGREAD_BINNING_PREPARATION.out.grouped_mappings
    )

    ch_multiqc_files = ch_multiqc_files.mix(SHORTREAD_BINNING_PREPARATION.out.bowtie2_assembly_multiqc)

    emit:
    grouped_mappings = ch_grouped_mappings
    versions         = ch_versions
    multiqc_files    = ch_multiqc_files
}
