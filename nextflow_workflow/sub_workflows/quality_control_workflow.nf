#!/usr/bin/env nextflow

include { qualityControl } from '../modules/qualityControl.nf'
include { qualityControl as qualityControlRepair } from '../modules/qualityControl.nf'
include { detectEncoding } from '../modules/detect_encoding.nf'
include { repairReads } from '../modules/repair_reads.nf'

workflow quality_control_workflow {
    take:
    quality_control_input_channel

    main:
    qualityControl(quality_control_input_channel)

    conditional_channel = qualityControl.out.qc_status_report.branch { _sample_id, exit_code, _fastqc_report ->
        success: exit_code == "0"
        failed: exit_code != "0"
    }
    quality_control_success_out_ch = qualityControl.out.quality_control_output.join(
        conditional_channel.success.map { sample_id, _exit_code, _fastqc_report -> sample_id }
    )

    input_read_channel = quality_control_input_channel.map { sample_id, read1, read2, _lower, _upper ->
        tuple(sample_id, read1, read2)
    }
    cutoffs_ch = quality_control_input_channel.map { sample_id, _read1, _read2, lower, upper ->
        tuple(sample_id, lower, upper)
    }

    // conditional_channel.failed.view { "Conditional Channel: ${it}" }
    // conditional_channel.failed
    conditional_channel.failed
        .map { sample_id, _exit_code, fastqc_report -> tuple(sample_id, fastqc_report) }
        .set { encoding_detection_input_ch }

    detectEncoding(encoding_detection_input_ch)
    detectEncoding.out.encoding_output_channel
    // .view { "Detected Encoding: ${it}" }

    input_read_channel
        .join(
            detectEncoding.out.encoding_output_channel
        )
        .set { repair_reads_input_ch }

    // repair_reads_input_ch.view { "Repair Reads Input: ${it}" }

    repairReads(repair_reads_input_ch)
    repairReads.out.repaired_reads_output_channel
    // .view { "Repaired Reads: ${it}" }

    repairReads.out.repaired_reads_output_channel
        .join(cutoffs_ch)
        .set { repaired_quality_control_input_ch }

    qualityControlRepair(repaired_quality_control_input_ch)
    quality_control_repair_out_ch = qualityControlRepair.out.quality_control_output
    quality_control_out_ch = quality_control_success_out_ch.mix(quality_control_repair_out_ch)

    emit:
    quality_control_out_ch
}
