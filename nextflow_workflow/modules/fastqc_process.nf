#!/usr/bin/env nextflow

process fastqc_process {

    tag "${sample_id}"

    publishDir "${parent_output_dir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(input_file1), path(input_file2)
    val parent_output_dir

    output:
    path "fastqc_output/*"

    script:
    """
    mkdir fastqc_output
    fastqc -o fastqc_output ${input_file1} ${input_file2}
    """
}
