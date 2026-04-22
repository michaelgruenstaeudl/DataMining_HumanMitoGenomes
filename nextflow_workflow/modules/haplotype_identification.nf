#!/usr/bin/env nextflow

process haplotype_identification {

    tag "${sample_id}"

    publishDir "${parent_output_dir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta_file)
    val phase
    val parent_output_dir

    output:
    tuple val(sample_id), path("${output_directory}/haplotype_output.csv"), emit: haplotype_output_channel
    path "${output_directory}/*"

    script:
    output_directory = "haplotype_identification/${phase}"

    """
    mkdir -p ${output_directory}
    haplogrep3 classify --input ${fasta_file} --out ${output_directory}/haplotype_output.csv  --tree phylotree-rcrs@17.2
    """
}
