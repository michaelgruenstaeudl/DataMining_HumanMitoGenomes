#!/usr/bin/env nextflow

process generate_sam {

    tag "${sample_id}"

    // container 'evolbioinfo/schmutzi:v1.5.6'
    publishDir "results/${sample_id}", mode: 'copy'
    container "quay.io/biocontainers/bowtie2:2.5.4--he96a11b_6"

    input:
    tuple val(sample_id), path(fastq1), path(fastq2), path(fasta_file)
    val phase

    output:
    tuple val(sample_id), path("sam_output_${phase}/${sample_id}.sam"), emit: bam_output_ch

    script:
    def output_sam_filename = "sam_output_${phase}/${sample_id}.sam"
    // def output_bam_filename = "bam_output/${sample_id}.bam"
    // #samtools view -@ 8 -Sb -T ${fasta_file} ${output_sam_filename} >${output_bam_filename}

    """
    ## MAPPING ORIGINAL READS
    mkdir -p sam_output_${phase}
    
    mkdir -p temp_data2/temp/db
    bowtie2-build ${fasta_file} temp_data2/temp/db/${sample_id}
    bowtie2 -x temp_data2/temp/db/${sample_id} -1 ${fastq1} -2 ${fastq2} -S ${output_sam_filename} -p 8
    """
}
