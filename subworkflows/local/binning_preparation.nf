
include { SHORTREAD_BINNING_PREPARATION  } from './shortread_binning_preparation'
include { LONGREAD_BINNING_PREPARATION    } from './longread_binning_preparation'

workflow BINNING_PREPARATION {
    take:
    shortread_assemblies           // channel: [ val(meta), path(assembly) ]
    shortreads                // channel: [ val(meta), [ reads ] ]
    longread_assemblies           // channel: [ val(meta), path(assembly) ]
    longreads                // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()
        // multiple symlinks to the same assembly -> use first of sorted list
    SHORTREAD_BINNING_PREPARATION ( shortread_assemblies, shortreads )
    LONGREAD_BINNING_PREPARATION ( longread_assemblies, longreads )

    ch_grouped_mappings = SHORTREAD_BINNING_PREPARATION.out.grouped_mappings
        .mix( LONGREAD_BINNING_PREPARATION.out.grouped_mappings )

    ch_versions = ch_versions.mix( SHORTREAD_BINNING_PREPARATION.out.versions )
    ch_versions = ch_versions.mix( LONGREAD_BINNING_PREPARATION.out.versions )
    ch_multiqc_files = ch_multiqc_files.mix( SHORTREAD_BINNING_PREPARATION.out.bowtie2_assembly_multiqc )

    emit:
    grouped_mappings         = ch_grouped_mappings
    versions                 = ch_versions
    multiqc_files            = ch_multiqc_files
}
