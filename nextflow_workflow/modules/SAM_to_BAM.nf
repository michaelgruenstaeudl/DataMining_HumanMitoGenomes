#!/usr/bin/env nextflow


process SAM_to_BAM {

    tag "${sample_id}"

    publishDir "results/${sample_id}", mode: 'copy', overwrite: true
    container "biocontainers/samtools:v1.9-4-deb_cv1"

    input:
    tuple val(sample_id), path(fasta_file), path(sam_file)
    val phase

    output:
    tuple val(sample_id), path("bam_output_${phase}/${sample_id}.bam"), emit: bam_output_ch

    script:
    def output_bam_filename = "bam_output_${phase}/${sample_id}.bam"

    """
    ## BAM FILE FOR ORIGINAL READS
    mkdir -p bam_output_${phase}
    samtools view -@ 8 -Sb -T ${fasta_file} ${sam_file} >${output_bam_filename}
    """
}
