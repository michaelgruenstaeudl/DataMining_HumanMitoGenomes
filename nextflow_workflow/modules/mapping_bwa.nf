#!/usr/bin/env nextflow

process mapping_bwa {

    label "bwa_env"

    tag "${sample_id}"

    input:
    // need to change publishDir to dsl=2 implementation, since dsl=2 does not support publishDir in process input block
    // publishDir "${parent_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    tuple val(sample_id), path(input_file1), path(input_file2), path(reference_fasta), path(index_file1), path(index_file2), path(index_file3), path(index_file4), path(index_file5)
    val parent_output_dir

    output:
    tuple val(sample_id), path("mapped_output/${sample_id}.sam"), emit: mapping_bwa_output

    script:
    """
    mkdir mapped_output
    bwa-mem2 mem -t 4 ${reference_fasta} ${input_file1} ${input_file2} > mapped_output/${sample_id}.sam
    """
}
