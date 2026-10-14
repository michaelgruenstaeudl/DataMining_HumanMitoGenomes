#!/usr/bin/env nextflow

process repair_reads {

    tag "${meta.id}"

    input:

    tuple val(meta), path(input_file_list), val(encoding_int)

    output:
    // val sample_id, emit: sample_id
    tuple val(meta), path("repair_read_output/${meta.id}_*.fastq", optional: true), emit: repaired_reads_output_channel
    tuple val(meta), path("repair_read_output"), optional: true, emit: repair_output_publish_dir

    script:
    def output_dir = "repair_read_output"
    def input_file = meta.single_end
        ? "in1=${input_file_list[0]}"
        : "in1=${input_file_list[0]} in2=${input_file_list[1]}"

    def output = meta.single_end
        ? "out=${output_dir}/${meta.id}_1.fastq outs=${output_dir}/${meta.id}_singletons.fastq"
        : "out=${output_dir}/${meta.id}_1.fastq out2=${output_dir}/${meta.id}_2.fastq outs=${output_dir}/${meta.id}_singletons.fastq"
    """
    mkdir ${output_dir}
    repair.sh \\
        ${input_file} \\
        ${output} \\
        overwrite=true \\
        qin=${encoding_int}
    """
}
