#!/usr/bin/env nextflow

process quality_control {

    tag "${sample_id}"

    // container 'community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18'
    errorStrategy 'ignore'
    publishDir "${parent_output_dir}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(input_file_list), val(is_single_end), val(lower_cutoff), val(upper_cutoff)
    val quality_control_stage
    val parent_output_dir

    output:
    // val sample_id, emit: sample_id
    tuple val(sample_id), env("exit_code"), path("${quality_control_stage}/${sample_id}_1.fastq_trimming_report.txt", optional: true), emit: qc_status_report
    tuple val(sample_id), path("${quality_control_stage}/${sample_id}_*.fq", optional: true), val(is_single_end), emit: quality_control_output
    path "${quality_control_stage}/*"

    script:
    // clip_length = upper_cutoff.toInteger() - 3
    //     trim_galore --paired_end --clip_R1 3 --clip_R2 3  --max_length ${upper_cutoff} --quality 20 \
    def input_file = is_single_end
        ? "${input_file_list[0]}"
        : "${input_file_list[0]} ${input_file_list[1]}"
    def paired_end_option = is_single_end ? "" : "--paired"
    """
    set +e
    trim_galore ${paired_end_option} --max_length ${upper_cutoff} --quality 20 \
    --fastqc \
    --output_dir ${quality_control_stage} \
    ${input_file}
    
    exit_code=\$(echo \$? )
    echo \$exit_code > .exitcode
    echo \$exit_code
    """
}
