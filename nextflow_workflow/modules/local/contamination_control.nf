#!/usr/bin/env nextflow

process contamination_removal {

    tag "${meta.id}"

    input:

    // container "community.wave.seqera.io/library/kraken2:2.1.6--62621491c438309a"

    tuple val(meta), path(fastq_input_list), path(db_path)
    val confidence_threshold

    output:
    tuple val(meta), path("kraken2_output/${meta.id}.unclassified*.fastq"), emit: contamination_removal_fastq_output
    tuple val(meta), path("kraken2_output"), emit: contamination_removal_publish_dir

    script:
    def paired_end_option = meta.single_end ? "" : "--paired"
    fastq_input_param = meta.single_end ? fastq_input_list[0] : "${fastq_input_list[0]} ${fastq_input_list[1]}"
    def unclassified_output = meta.single_end ? "${meta.id}.unclassified.fastq" : "${meta.id}.unclassified#.fastq"
    """
    mkdir kraken2_output
    kraken2 \\
        --db ${db_path} \\
        ${paired_end_option} \\
        --confidence ${confidence_threshold}\\
        --report kraken2_output/${meta.id}_kraken2_report.txt  \\
        --unclassified-out kraken2_output/"${unclassified_output}" \\
        --output kraken2_output/${meta.id}_kraken2_output.kraken \\
        ${fastq_input_param}
    """
}
