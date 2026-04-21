#!/usr/bin/env nextflow

process determine_haplotype {

    tag "${sample_id}"

    publishDir "${parent_output_dir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta_file)
    val parent_output_dir
    val phase

    output:
    tuple val(sample_id), path("${phase}/haplotype_output.csv"), emit: haplogrep_output_channel

    script:
    """
    mkdir ${phase}
    haplogrep3 classify --input ${fasta_file} --out ${phase}/haplotype_output.csv  --tree phylotree-rcrs@17.2
    """
}
