nextflow_workflow {
    name "Test Subworkflow GENERATE_DOWNSTREAM_SAMPLESHEETS"
    script "../main.nf"
    workflow "GENERATE_DOWNSTREAM_SAMPLESHEETS"

    tag "subworkflows"
    tag "subworkflows_local"
    tag "subworkflows/generate_downstream_samplesheets"

    test("reads - taxprofiler,funscan") {

        when {
            params {
                modules_testdata_base_path      = "https://raw.githubusercontent.com/nf-core/test-datasets/"
                outdir                          = "."
                generate_pipeline_samplesheets  = 'taxprofiler,funscan'
            }
            workflow {
                """
                input[0] = Channel.of(
                        [
                            [id:'test_taxprofiler_funscan', single_end:false, long_reads:true, amount_of_files:3],
                            file(params.modules_testdata_base_path + 'mag/samplesheets/samplesheet.hybrid.csv', checkIfExists: true)
                        ]
                )
                """
            }
        }

        then {
            assertAll(
                { assert workflow.success},
                { assert snapshot(
                    [
                        "${params.outdir}/downstream_samplesheets/funscan.csv",
                        "${params.outdir}/downstream_samplesheets/taxprofiler.csv"
                    ]).match()
                },
            )
        }
    }
}
