#!/usr/bin/env nextflow

process fastqc_process {

    tag "${meta.id}"

    input:

    tuple val(meta), path(input_file_list)
    val phase

    output:
    tuple val(meta), path("${phase}/${input_file_list[0].simpleName}_fastqc.zip"), emit: fastqc_output_channel
    tuple val(meta), path("${phase}"), emit: fastqc_publish_dir

    script:
    def input_file1 = input_file_list[0]
    def input_file2 = meta.single_end ? null : input_file_list[1]
    """
    set -euo pipefail
    mkdir ${phase}
    fastqc -o ${phase} ${input_file1} ${input_file2}
    """
}
