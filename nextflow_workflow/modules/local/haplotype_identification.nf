#!/usr/bin/env nextflow

process haplotype_identification {

    tag "${meta.id}"

    input:
    tuple val(meta), path(fasta_file)
    val phase

    output:
    tuple val(meta), path("${output_directory}/haplotype_output.csv"), emit: haplotype_output_channel

    script:
    // tuple val(meta), path("${output_directory}"), emit: haplotype_output_publish_dir
    output_directory = "haplotype_identification/${phase}"

    """
    set -euo pipefail
    mkdir -p ${output_directory}
    haplogrep3 classify --input ${fasta_file} --out ${output_directory}/haplotype_output.csv  --tree phylotree-rcrs@17.2
    """
}
