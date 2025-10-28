#!/usr/bin/env nextflow

include { qualityControl } from '../modules/qualityControl.nf'
include { qualityControl as qualityControlRepair } from '../modules/qualityControl.nf'
include { detectEncoding } from '../modules/detect_encoding.nf'
include { repairReads } from '../modules/repair_reads.nf'
include { addSizeTracking ; WriteCSV } from '../lib/utils.nf'
include { calculate_sequence_length_threshold } from '../modules/compute_sequence_length_threshold.nf'
include { calculate_sequence_length_threshold as calculate_sequence_length_threshold_repair } from '../modules/compute_sequence_length_threshold.nf'
workflow quality_control_workflow {
    take:
    quality_control_input_channel
    calculate_sequence_length_threshold_script_ch

    main:
    repair_read_size_ch = channel.empty()

    // Calculate sequence length threshold process
    input_to_calculate_sequence_length_threshold = quality_control_input_channel
        .map { sample_id, read1, _read2 ->
            tuple(sample_id, read1)
        }
        .combine(
            calculate_sequence_length_threshold_script_ch
        )
    calculate_sequence_length_threshold(input_to_calculate_sequence_length_threshold)

    // Initial Quality control process
    quality_control_input_channel_with_cutoffs = quality_control_input_channel.join(calculate_sequence_length_threshold.out.length_cutoffs)
    qualityControl(quality_control_input_channel_with_cutoffs, 'initial_quality_control')

    // Separate channels for success and failed QC based on exit code
    conditional_channel = qualityControl.out.qc_status_report.branch { _sample_id, exit_code, _fastqc_report ->
        success: exit_code == "0"
        failed: exit_code != "0"
    }
    quality_control_success_out_ch = qualityControl.out.quality_control_output.join(
        conditional_channel.success.map { sample_id, _exit_code, _fastqc_report -> sample_id }
    )
    success_out_cutoff_ch = calculate_sequence_length_threshold.out.length_cutoffs.join(
        conditional_channel.success.map { sample_id, _exit_code, _fastqc_report -> sample_id }
    )

    // Encoding detection process for failed QC samples
    conditional_channel.failed
        .map { sample_id, _exit_code, fastqc_report -> tuple(sample_id, fastqc_report) }
        .set { encoding_detection_input_ch }

    detectEncoding(encoding_detection_input_ch)


    // Repair reads process for failed QC samples
    quality_control_input_channel
        .join(
            detectEncoding.out.encoding_output_channel
        )
        .set { repair_reads_input_ch }

    repairReads(repair_reads_input_ch, "qc")

    repair_read_size_ch = addSizeTracking(repair_read_size_ch, repairReads.out.repaired_reads_output_channel, "Repair Read Output")

    // Calculate sequence length threshold for repaired reads
    input_to_calculate_sequence_length_threshold_repair = repairReads.out.repaired_reads_output_channel
        .map { sample_id, read1, _read2 ->
            tuple(sample_id, read1)
        }
        .combine(
            calculate_sequence_length_threshold_script_ch
        )
    calculate_sequence_length_threshold_repair(input_to_calculate_sequence_length_threshold_repair)

    // Quality control process for repaired reads
    repairReads.out.repaired_reads_output_channel
        .join(calculate_sequence_length_threshold_repair.out.length_cutoffs)
        .set { repaired_quality_control_input_ch }

    qualityControlRepair(repaired_quality_control_input_ch, 'repaired_quality_control')
    quality_control_repair_out_ch = qualityControlRepair.out.quality_control_output

    // Merging quality control outputs and cutoff outputs from both initial(success) and repaired reads
    quality_control_out_ch = quality_control_success_out_ch.mix(quality_control_repair_out_ch)

    cutoffs_ch = calculate_sequence_length_threshold_repair.out.length_cutoffs.mix(success_out_cutoff_ch)

    emit:
    quality_control_out_ch
    repair_read_size_ch
    cutoffs_ch
}
