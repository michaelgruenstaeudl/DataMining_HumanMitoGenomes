#!/usr/bin/env nextflow

process fastqc_process {

    tag "${sample_id}"

    publishDir "${parent_output_dir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(input_file_list), val(single_end)
    val parent_output_dir
    val phase

    output:
    path "${phase}/*"

    script:
    def input_file1 = input_file_list[0]
    def input_file2 = single_end ? null : input_file_list[1]
    """
    mkdir ${phase}
    fastqc -o ${phase} ${input_file1} ${input_file2}
    """
}
