#!/usr/bin/env nextflow

process contamination_removal {

    tag "${sample_id}"

    // container 'evolbioinfo/schmutzi:v1.5.6'
    publishDir "results/${sample_id}", mode: 'copy'
    container "community.wave.seqera.io/library/kraken2:2.1.6--62621491c438309a"

    input:
    tuple val(sample_id), path(db_path), path(fastq1), path(fastq2)

    output:
    tuple val(sample_id), path("kraken2_output/${sample_id}.unclassified_1.fastq"), path("kraken2_output/${sample_id}.unclassified_2.fastq"), emit: contamination_removal_fastq_output
    path "kraken2_output/${sample_id}_kraken2_report.txt"
    path "kraken2_output/${sample_id}_kraken2_output.kraken"

    script:
    // """
    // # schmutzi.pl --ref ${reference_data} --t 8 results /projects1/tools/schmutzi/alleleFreqMT/197/freqs/ ${input_bam_file}

    // """
    """
    mkdir kraken2_output
    kraken2 \\
        --db /home/b_thapamagar/BioInformatics/DataMining_HumanMitoGenomes/nextflow_workflow/db \\
        --paired SRR6245223_1.fastq SRR6245223_2.fastq \\
        --report kraken2_output/${sample_id}_kraken2_report.txt  \\
        --unclassified-out kraken2_output/"${sample_id}."unclassified#.fastq \\
        --output kraken2_output/${sample_id}_kraken2_output.kraken
    """
}
