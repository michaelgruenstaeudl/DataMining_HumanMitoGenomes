#!/usr/bin/env nextflow

process quality_control {

    tag "${meta.id}"

    // container 'community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18'
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(input_file_list), val(lower_cutoff), val(upper_cutoff)
    val quality_control_stage

    output:
    // val sample_id, emit: sample_id
    tuple val(meta), env("exit_code"), path("${quality_control_stage}/${meta.id}*.fastq_trimming_report.txt", optional: true), emit: qc_status_report
    tuple val(meta), path("${quality_control_stage}/${meta.id}_*.fq", optional: true), emit: quality_control_output
    tuple val(meta), path("${quality_control_stage}"), emit: quality_control_publish_dir

    script:
    // clip_length = upper_cutoff.toInteger() - 3
    //     trim_galore --paired_end --clip_R1 3 --clip_R2 3  --max_length ${upper_cutoff} --quality 20 \
    def input_file = meta.single_end
        ? "${input_file_list[0]}"
        : "${input_file_list[0]} ${input_file_list[1]}"
    def paired_end_option = meta.single_end ? "" : "--paired"
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
