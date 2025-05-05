#!/usr/bin/env nextflow

process mapping_process {

    tag "${sample_id}"

    publishDir 'results', mode: 'copy'
    container 'fastq_sifter'

    input:
    tuple val(sample_id), path(input_file1), path(input_file2)
    path reference_fasta

    output:
    path "mapped_output/*"

    script:
    """
    mkdir mapped_output
    FastqSifter --out mapped_output/map --FASTA=${reference_fasta} --left=${input_file1} --right=${input_file2}
    """
}
