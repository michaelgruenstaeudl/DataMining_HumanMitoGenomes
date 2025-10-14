#!/usr/bin/env nextflow

process qualityControl {

    tag "${sample_id}"
    errorStrategy 'ignore'
    publishDir "results/${sample_id}", mode: 'copy', overwrite: true
    container 'community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18'

    input:
    tuple val(sample_id), path(input_file1), path(input_file2), val(lower_cutoff), val(upper_cutoff)
    val quality_control_stage

    output:
    // val sample_id, emit: sample_id
    tuple val(sample_id), env("exit_code"), path("${quality_control_stage}/${sample_id}_1.fastq_trimming_report.txt", optional: true), emit: qc_status_report
    tuple val(sample_id), path("${quality_control_stage}/${sample_id}_1_val_1.fq", optional: true), path("${quality_control_stage}/${sample_id}_2_val_2.fq", optional: true), emit: quality_control_output
    path "${quality_control_stage}/${sample_id}_1_val_1_fastqc.html", optional: true
    path "${quality_control_stage}/${sample_id}_1_val_1_fastqc.zip", optional: true
    path "${quality_control_stage}/${sample_id}_1.fastq_trimming_report.txt", optional: true, emit: trimming_report_1
    path "${quality_control_stage}/${sample_id}_2_val_2_fastqc.html", optional: true
    path "${quality_control_stage}/${sample_id}_2_val_2_fastqc.zip", optional: true
    path "${quality_control_stage}/${sample_id}_2.fastq_trimming_report.txt", optional: true

    script:
    clip_length = upper_cutoff.toInteger() - 3
    """
    set +e
    trim_galore --paired_end --three_prime_clip_R1 3 --three_prime_clip_R2 3 --length ${clip_length} --max_length ${upper_cutoff} --quality 10 \
    --fastqc \
    --output_dir ${quality_control_stage} \
    ${input_file1} ${input_file2} 
    
    exit_code=\$(echo \$? )
    echo \$exit_code > .exitcode
    echo \$exit_code
    """
}
