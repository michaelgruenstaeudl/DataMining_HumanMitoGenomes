#!/usr/bin/env nextflow


process sam_to_bam {

    tag "${sample_id}"

    // container "biocontainers/samtools:v1.9-4-deb_cv1" --docker container
    // container 'oras://community.wave.seqera.io/library/samtools:1.22.1--9a10f06c24cdf05f' --apptainer container

    publishDir "${parent_output_dir}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(fasta_file), path(sam_file)
    val phase
    val parent_output_dir

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
