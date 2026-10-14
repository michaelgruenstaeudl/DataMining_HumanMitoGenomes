#!/usr/bin/env nextflow

process fastqc_process {

    tag "${sample_id}"

    input:

    tuple val(sample_id), path(input_file_list), val(single_end)
    val phase

    output:
    tuple val(sample_id), path("${phase}/${input_file_list[0].simpleName}_fastqc.zip"), emit: fastqc_output_channel
    tuple val(sample_id), path("${phase}"), emit: fastqc_publish_dir

    script:
    def input_file1 = input_file_list[0]
    def input_file2 = single_end ? null : input_file_list[1]
    """
    mkdir ${phase}
    fastqc -o ${phase} ${input_file1} ${input_file2}
    """
}
