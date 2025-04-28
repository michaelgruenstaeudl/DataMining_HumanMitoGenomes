#!/usr/bin/env nextflow

process qualityControl {

    tag "${sample_id}"

    publishDir 'results', mode: 'copy'
    container 'community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18'

    input:
    tuple val(sample_id), path(input_file1), path(input_file2)
    tuple val(lower_cutoff), val(upper_cutoff)

    output:
    path "trimmed_output/${sample_id}_1_val_1_fastqc.html", emit: fastqc_report1
    path "trimmed_output/${sample_id}_1_val_1_fastqc.zip", emit: fastqc_report1_zip
    path "trimmed_output/${sample_id}_1_val_1.fq", emit: trimmed_fastq1
    path "trimmed_output/${sample_id}_1.fastq_trimming_report.txt", emit: fastq_trimming_report1
    path "trimmed_output/${sample_id}_2_val_2_fastqc.html", emit: fastqc_report2
    path "trimmed_output/${sample_id}_2_val_2_fastqc.zip", emit: fastqc_report2_zip
    path "trimmed_output/${sample_id}_2_val_2.fq", emit: trimmed_fastq2
    path "trimmed_output/${sample_id}_2.fastq_trimming_report.txt", emit: fastq_trimming_report2

    script:
    """
    trim_galore --paired_end --three_prime_clip_R1 3 --three_prime_clip_R2 3 --length 72 --max_length ${upper_cutoff} --quality 20 \
    --fastqc \
    --output_dir trimmed_output \
    ${input_file1} ${input_file2}
    """
}
