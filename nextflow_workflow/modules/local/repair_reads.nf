#!/usr/bin/env nextflow

process repair_reads {

    tag "${sample_id}"

    input:

    tuple val(sample_id), path(input_file_list), val(is_single_end), val(encoding_int)

    output:
    // val sample_id, emit: sample_id
    tuple val(sample_id), path("repair_read_output/${sample_id}_*.fastq", optional: true), val(is_single_end), emit: repaired_reads_output_channel
    tuple val(sample_id), path("repair_read_output"), optional: true, emit: repair_output_publish_dir

    script:
    def output_dir = "repair_read_output"
    def input_file = is_single_end
        ? "in1=${input_file_list[0]}"
        : "in1=${input_file_list[0]} in2=${input_file_list[1]}"

    def output = is_single_end
        ? "out=${output_dir}/${sample_id}_1.fastq outs=${output_dir}/${sample_id}_singletons.fastq"
        : "out=${output_dir}/${sample_id}_1.fastq out2=${output_dir}/${sample_id}_2.fastq outs=${output_dir}/${sample_id}_singletons.fastq"
    """
    mkdir ${output_dir}
    repair.sh \\
        ${input_file} \\
        ${output} \\
        overwrite=true \\
        qin=${encoding_int}
    """
}
