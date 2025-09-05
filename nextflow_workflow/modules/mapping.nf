#!/usr/bin/env nextflow

process mapping_process {

    tag "${sample_id}"

    publishDir "results/${sample_id}", mode: 'copy'
    container 'fastq_sifter'

    input:
    tuple val(sample_id), path(input_file1), path(input_file2), path(reference_fasta)

    output:
    tuple val(sample_id), path("mapped_output/${sample_id}.filtered.A.fq"), path("mapped_output/${sample_id}.filtered.B.fq"), path("mapped_output/${sample_id}.bwa.A.sai"), path("mapped_output/${sample_id}.bwa.B.sai"), path("mapped_output/${sample_id}.bwa.sampe.sam"), emit: mapping_process_output

    script:
    """
    mkdir mapped_output
    FastqSifter --out mapped_output/${sample_id} --FASTA=${reference_fasta} --left=${input_file1} --right=${input_file2}
    """
}
