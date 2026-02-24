#!/usr/bin/env nextflow

process contamination_removal {

    tag "${sample_id}"

    // container "community.wave.seqera.io/library/kraken2:2.1.6--62621491c438309a"

    publishDir "${parent_output_dir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq_input_list), val(is_single_end), path(db_path)
    val confidence_threshold
    val parent_output_dir

    output:
    tuple val(sample_id), path("kraken2_output/${sample_id}.unclassified*.fastq"), val(is_single_end), emit: contamination_removal_fastq_output
    path "kraken2_output/${sample_id}_kraken2_report.txt"
    path "kraken2_output/${sample_id}_kraken2_output.kraken"

    script:
    def paired_end_option = is_single_end ? "" : "--paired"
    fastq_input_param = is_single_end ? fastq_input_list[0] : "${fastq_input_list[0]} ${fastq_input_list[1]}"
    def unclassified_output = is_single_end ? "${sample_id}.unclassified.fastq" : "${sample_id}.unclassified#.fastq"
    """
    mkdir kraken2_output
    kraken2 \\
        --db ${db_path} \\
        ${paired_end_option} \\
        --confidence ${confidence_threshold}\\
        --report kraken2_output/${sample_id}_kraken2_report.txt  \\
        --unclassified-out kraken2_output/"${unclassified_output}" \\
        --output kraken2_output/${sample_id}_kraken2_output.kraken \\
        ${fastq_input_param}
    """
}
