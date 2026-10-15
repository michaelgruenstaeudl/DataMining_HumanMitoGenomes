#!/usr/bin/env nextflow


process sam_to_bam {

    label "samtools_env"
    tag "${meta.id}"

    input:
    tuple val(meta), path(fasta_file), path(sam_file)
    val phase

    output:
    tuple val(meta), path("bam_output_${phase}/${meta.id}.bam"), emit: bam_output_ch
    tuple val(meta), path("bam_output_${phase}"), emit: bam_output_publish_dir

    script:
    def output_bam_filename = "bam_output_${phase}/${meta.id}.bam"

    """
    set -euo pipefail
    ## BAM FILE FOR ORIGINAL READS
    mkdir -p bam_output_${phase}
    samtools view -@ 8 -Sb -T ${fasta_file} ${sam_file} >${output_bam_filename}
    """
}
