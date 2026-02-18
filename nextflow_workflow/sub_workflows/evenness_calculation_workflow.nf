#!/usr/bin/env nextflow

include { generate_sam } from '../modules/generate_sam.nf'
include { sam_to_bam } from '../modules/sam_to_bam.nf'
include { calculate_evenness } from '../modules/calculate_evenness.nf'
include { repair_reads } from '../modules/repair_reads.nf'

workflow evenness_calculation_workflow {
    take:
    evenness_calculation_input_channel
    phase
    parent_output_dir

    main:
    fasta_file_ch = evenness_calculation_input_channel.map { sample_id, _fastq1, _fastq2, _gb_file, fasta_file ->
        tuple(sample_id, fasta_file)
    }


    gb_file_ch = evenness_calculation_input_channel.map { sample_id, _fastq1, _fastq2, gb_file, _fasta_file ->
        tuple(sample_id, gb_file)
    }
    generate_sam_input_ch = channel.empty()
    if (phase == 'initial_phase') {
        repair_reads_input_ch = evenness_calculation_input_channel
            .map { sample_id, fastq1, fastq2, _gb_file, _fasta_file ->
                tuple(sample_id, fastq1, fastq2)
            }
            .combine(channel.value("33"))

        repair_reads(repair_reads_input_ch, phase, parent_output_dir)
        repair_reads.out.repaired_reads_output_channel.join(fasta_file_ch).set { generate_sam_input_ch }
    }
    else {
        generate_sam_input_ch = evenness_calculation_input_channel
            .map { sample_id, fastq1, fastq2, _gb_file, _fasta_file ->
                tuple(sample_id, fastq1, fastq2)
            }
            .join(fasta_file_ch)
    }

    generate_sam(generate_sam_input_ch, phase, parent_output_dir)

    sam_to_bam_input_ch = fasta_file_ch.join(generate_sam.out.sam_output_ch)
    sam_to_bam(sam_to_bam_input_ch, phase, parent_output_dir)
    calculate_evenness_input_ch = gb_file_ch.join(sam_to_bam.out.bam_output_ch)
    calculate_evenness(calculate_evenness_input_ch, phase, parent_output_dir)
    evenness_output_ch = calculate_evenness.out.evenness_output_ch

    emit:
    evenness_output_ch
}
