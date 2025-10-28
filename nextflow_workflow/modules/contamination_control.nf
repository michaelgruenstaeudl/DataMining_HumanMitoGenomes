#!/usr/bin/env nextflow

process contamination_removal {

    tag "${sample_id}"

    // container "community.wave.seqera.io/library/kraken2:2.1.6--62621491c438309a"

    publishDir "results/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(fastq1), path(fastq2), path(db_path)

    output:
    tuple val(sample_id), path("kraken2_output/${sample_id}.unclassified_1.fastq"), path("kraken2_output/${sample_id}.unclassified_2.fastq"), emit: contamination_removal_fastq_output
    path "kraken2_output/${sample_id}_kraken2_report.txt"
    path "kraken2_output/${sample_id}_kraken2_output.kraken"

    script:
    """
    mkdir kraken2_output
    kraken2 \\
        --db ${db_path} \\
        --paired ${fastq1} ${fastq2} \\
        --report kraken2_output/${sample_id}_kraken2_report.txt  \\
        --unclassified-out kraken2_output/"${sample_id}."unclassified#.fastq \\
        --output kraken2_output/${sample_id}_kraken2_output.kraken
    """
}
