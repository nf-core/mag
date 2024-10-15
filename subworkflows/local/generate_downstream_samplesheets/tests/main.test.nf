nextflow_workflow {
    name "Test Subworkflow GENERATE_DOWNSTREAM_SAMPLESHEETS"
    script "../main.nf"
    workflow "GENERATE_DOWNSTREAM_SAMPLESHEETS"

    tag "subworkflows"
    tag "subworkflows_local"
    tag "subworkflows/generate_downstream_samplesheets"

    test("reads,assemblies - taxprofiler,funcscan") {

        when {
            params {
                modules_testdata_base_path      = "https://github.com/nf-core/test-datasets/raw/mag/test_data/"
                outdir                          = "."
                generate_pipeline_samplesheets  = 'taxprofiler,funcscan'
            }
            workflow {
                """
                input[0] = Channel.of(
                        [
                            [id:'test_minigut', group:0, single_end:false, amount_of_files:2],
                            file(params.modules_testdata_base_path + 'mag/test_data/test_minigut_R1.fastq.gz', checkIfExists: true)
                            file(params.modules_testdata_base_path + 'mag/test_data/test_minigut_R2.fastq.gz', checkIfExists: true)
                        ],
                        [
                            [id:'test_minigut_sample2', group:0, single_end:false, amount_of_files:2],
                            file(params.modules_testdata_base_path + 'mag/test_data/test_minigut_sample2_R1.fastq.gz', checkIfExists: true)
                            file(params.modules_testdata_base_path + 'mag/test_data/test_minigut_sample2_R2.fastq.gz', checkIfExists: true)
                        ]
                )
                input[1] = Channel.of(
                        [
                            [id:'test_minigut_spades', group:0, single_end:false, assembler:SPAdes, amount_of_files:1],
                            file(params.modules_testdata_base_path + 'mag/assemblies/SPAdes-test_minigut_contigs.fasta.gz', checkIfExists: true)
                        ],
                        [
                            [id:'test_minigut_megahit', group:0, single_end:false, assembler:MEGAHIT, amount_of_files:1],
                            file(params.modules_testdata_base_path + 'mag/assemblies/MEGAHIT-test_minigut.contigs.fa.gz', checkIfExists: true)
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
                        "${params.outdir}/downstream_samplesheets/funcscan.csv",
                        "${params.outdir}/downstream_samplesheets/taxprofiler.csv"
                    ]).match()
                },
            )
        }
    }
}
