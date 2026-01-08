#!/usr/bin/env nextflow

include { mapping_bwa } from '../modules/mapping_bwa.nf'
include { extract_mapped_fastq } from '../modules/extract_mapped_fastq.nf'

workflow mapping_workflow {
    take:
    mapping_input_ch
    parent_output_dir

    main:
    // generate_sam_input_ch = evenness_calculation_input_channel.map { sample_id, fastq1, fastq2, _gb_file, fasta_file ->
    //     tuple(sample_id, fastq1, fastq2, fasta_file)
    // }
    mapping_input_ch.view { it -> "Mapping input: ${it}" }
    mapping_bwa(mapping_input_ch, parent_output_dir)
    extract_mapped_fastq(mapping_bwa.out.mapping_bwa_output, parent_output_dir)
    mapping_process_output = extract_mapped_fastq.out.extract_mapped_fastq_output_ch

    emit:
    mapping_process_output
}
