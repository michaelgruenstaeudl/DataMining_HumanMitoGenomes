#!/usr/bin/env nextflow

process repair_reads {

    tag "${sample_id}"

    // container 'community.wave.seqera.io/library/bbmap:39.33--60639c9e1473b7a8'
    // container 'oras://community.wave.seqera.io/library/bbmap:39.37--9bf150ff9855d6f6'

    publishDir "${parent_output_dir}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(input_file_list), val(is_single_end), val(encoding_int)
    val phase
    val parent_output_dir

    output:
    // val sample_id, emit: sample_id
    tuple val(sample_id), path("repair_read_output_${phase}/${sample_id}_*.fastq", optional: true), val(is_single_end), emit: repaired_reads_output_channel
    path "repair_read_output_${phase}/${sample_id}_singletons.fastq", optional: true

    script:
    def output_dir = "repair_read_output_${phase}"
    def input_file = is_single_end
        ? "in1=${input_file_list[0]}"
        : "in1=${input_file_list[0]} in2=${input_file_list[1]}"

    def output = is_single_end
        ? "out1=${output_dir}/${sample_id}_1.fastq outs=${output_dir}/${sample_id}_singletons.fastq"
        : "out1=${output_dir}/${sample_id}_1.fastq out2=${output_dir}/${sample_id}_2.fastq outs=${output_dir}/${sample_id}_singletons.fastq"
    """
    mkdir ${output_dir}
    repair.sh \\
        ${input_file} \\
        ${output} \\
        overwrite=true \\
        qin=${encoding_int}
    """
}
