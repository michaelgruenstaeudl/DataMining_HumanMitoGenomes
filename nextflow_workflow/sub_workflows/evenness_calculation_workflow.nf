#!/usr/bin/env nextflow

include { generate_sam } from '../modules/generate_sam.nf'
include { SAM_to_BAM } from '../modules/SAM_to_BAM.nf'
include { calculate_evenness } from '../modules/calculate_evenness.nf'


workflow evenness_calculation_workflow {
    take:
    evenness_calculation_input_channel
    phase

    main:
    generate_sam_input_ch = evenness_calculation_input_channel.map { sample_id, fastq1, fastq2, _gb_file, fasta_file ->
        tuple(sample_id, fastq1, fastq2, fasta_file)
    }
    fasta_file_ch = evenness_calculation_input_channel.map { sample_id, _fastq1, _fastq2, _gb_file, fasta_file ->
        tuple(sample_id, fasta_file)
    }


    gb_file_ch = evenness_calculation_input_channel.map { sample_id, _fastq1, _fastq2, gb_file, _fasta_file ->
        tuple(sample_id, gb_file)
    }

    generate_sam(generate_sam_input_ch, phase)
    // // generate_sam.out.bam_output_ch.view { "Generated SAM files: ${it}" }
    sam_to_bam_input_ch = fasta_file_ch.join(generate_sam.out.bam_output_ch)
    SAM_to_BAM(sam_to_bam_input_ch, phase)
    // SAM_to_BAM.out.bam_output_ch.view { "Generated BAM files: ${it}" }
    calculate_evenness_input_ch = gb_file_ch.join(SAM_to_BAM.out.bam_output_ch)
    calculate_evenness(calculate_evenness_input_ch, phase)
    evenness_output_ch = calculate_evenness.out.evenness_output_ch

    emit:
    evenness_output_ch
}
