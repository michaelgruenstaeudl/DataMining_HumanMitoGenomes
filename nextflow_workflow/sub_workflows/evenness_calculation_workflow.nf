#!/usr/bin/env nextflow

include { generate_sam } from '../modules/generate_sam.nf'
include { sam_to_bam } from '../modules/sam_to_bam.nf'
include { calculate_evenness } from '../modules/calculate_evenness.nf'
include { repairReads } from '../modules/repair_reads.nf'

workflow evenness_calculation_workflow {
    take:
    evenness_calculation_input_channel
    phase

    main:
    // generate_sam_input_ch = evenness_calculation_input_channel.map { sample_id, fastq1, fastq2, _gb_file, fasta_file ->
    //     tuple(sample_id, fastq1, fastq2, fasta_file)
    // }
    fasta_file_ch = evenness_calculation_input_channel.map { sample_id, _fastq1, _fastq2, _gb_file, fasta_file ->
        tuple(sample_id, fasta_file)
    }


    gb_file_ch = evenness_calculation_input_channel.map { sample_id, _fastq1, _fastq2, gb_file, _fasta_file ->
        tuple(sample_id, gb_file)
    }
    repair_reads_input_ch = evenness_calculation_input_channel
        .map { sample_id, fastq1, fastq2, _gb_file, _fasta_file ->
            tuple(sample_id, fastq1, fastq2)
        }
        .combine(channel.value("33"))

    //encoding integer hardcoded to 33 for now

    repairReads(repair_reads_input_ch, phase)
    repairReads.out.repaired_reads_output_channel.join(fasta_file_ch).set { generate_sam_input_ch }
    generate_sam(generate_sam_input_ch, phase)

    // Separate channels for success and failed generate_sam based on exit code
    //     conditional_channel = generate_sam.out.sam_status_report.branch { _sample_id, exit_code ->
    //         success: exit_code == "0"
    //         failed: exit_code != "0"
    //     }
    //     generate_sam_success_out_ch = generate_sam.out.sam_output_ch.join(
    //         conditional_channel.success.map { sample_id, _exit_code, _fastqc_report -> sample_id }
    //     )

    //     generate_sam_input_ch
    //         .join(conditional_channel.failed.map { sample_id, _exit_code -> tuple(sample_id) })
    //         .set { repair_input_ch }

    // generate_sam
    // generate_sam.out.bam_output_ch.view { "Generated SAM files: ${it}" }
    sam_to_bam_input_ch = fasta_file_ch.join(generate_sam.out.sam_output_ch)
    sam_to_bam(sam_to_bam_input_ch, phase)
    // SAM_to_BAM.out.bam_output_ch.view { "Generated BAM files: ${it}" }
    calculate_evenness_input_ch = gb_file_ch.join(sam_to_bam.out.bam_output_ch)
    calculate_evenness(calculate_evenness_input_ch, phase)
    evenness_output_ch = calculate_evenness.out.evenness_output_ch

    emit:
    evenness_output_ch
}
