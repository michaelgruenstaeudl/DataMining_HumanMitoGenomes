#!/usr/bin/env nextflow

include { quality_control } from '../modules/quality_control.nf'
// include { quality_control as quality_control_repair } from '../modules/quality_control.nf'
include { detectEncoding } from '../modules/detect_encoding.nf'
include { repair_reads } from '../modules/repair_reads.nf'
include { addSizeTracking } from '../lib/utils.nf'
include { calculate_sequence_length_threshold } from '../modules/compute_sequence_length_threshold.nf'
// include { calculate_sequence_length_threshold as calculate_sequence_length_threshold_repair } from '../modules/compute_sequence_length_threshold.nf'
include { fastqc_process ; fastqc_process as fastqc_process_after_qc } from '../modules/fastqc_process.nf'
// include { quality_control_cutadapt as quality_control ; quality_control_cutadapt as quality_control_repair } from '../modules/quality_control_cutadapt.nf'
workflow quality_control_workflow {
    take:
    quality_control_input_channel
    calculate_sequence_length_threshold_script_ch
    parent_output_dir

    main:
    repair_read_size_ch = channel.empty()

    // quality_control_input_channel.view { it -> "Quality control input: ${it}" }
    //initial fastqc report 
    fastqc_process(quality_control_input_channel, parent_output_dir, "initial_fastqc")
    initial_fastqc_output_ch = fastqc_process.out.fastqc_output_channel
    detectEncoding(initial_fastqc_output_ch)

    // Calculate sequence length threshold process
    // input_to_calculate_sequence_length_threshold = quality_control_input_channel
    //     .map { sample_id, reads, _is_single_end ->
    //         tuple(sample_id, reads[0])
    //     }
    //     .combine(
    //         calculate_sequence_length_threshold_script_ch
    //     )
    // calculate_sequence_length_threshold(input_to_calculate_sequence_length_threshold)

    // Initial Quality control process
    // quality_control_input_channel_with_cutoffs = quality_control_input_channel.join(calculate_sequence_length_threshold.out.length_cutoffs)
    // quality_control(quality_control_input_channel_with_cutoffs, 'initial_quality_control', parent_output_dir)
    // quality_control_output = quality_control.out.quality_control_output.map { sample_id, reads, is_single_end ->
    //     def reads_list = reads instanceof List
    //         ? reads.findAll { item ->
    //             item.name =~ /_(val_1|val_2)\.fq$/
    //         }
    //         : [reads]
    //     tuple(sample_id, reads_list, is_single_end)
    // }
    // // quality_control_output.view { it -> "Quality control output: ${it}" }
    // // Separate channels for success and failed QC based on exit code
    // qc_status_report = quality_control.out.qc_status_report.map { sample_id, exit_code, fastqc_report ->
    //     def fastqc_report_item = fastqc_report instanceof List ? fastqc_report[0] : fastqc_report
    //     tuple(sample_id, exit_code, fastqc_report_item)
    // }
    // conditional_channel = qc_status_report.branch { _sample_id, exit_code, _fastqc_report ->
    //     success: exit_code == "0"
    //     failed: exit_code != "0"
    // }
    // quality_control_success_out_ch = quality_control_output.join(
    //     conditional_channel.success.map { sample_id, _exit_code, _fastqc_report -> sample_id }
    // )
    // // quality_control_success_out_ch.view { it -> "Quality control success output: ${it}" }
    // success_out_cutoff_ch = calculate_sequence_length_threshold.out.length_cutoffs.join(
    //     conditional_channel.success.map { sample_id, _exit_code, _fastqc_report -> sample_id }
    // )

    // // Encoding detection process for failed QC samples
    // conditional_channel.failed
    //     .map { sample_id, _exit_code, fastqc_report -> tuple(sample_id, fastqc_report) }
    //     .set { encoding_detection_input_ch }

    // detectEncoding(encoding_detection_input_ch)
    // detectEncoding.out.encoding_output_channel.view { it -> "Encoding detection output: ${it}" }

    // Repair reads process for failed QC samples
    repair_reads_input_ch = quality_control_input_channel.join(
        detectEncoding.out.encoding_output_channel
    )

    repair_reads(repair_reads_input_ch, parent_output_dir)
    repair_read_output_ch = repair_reads.out.repaired_reads_output_channel.map { sample_id, reads, is_single_end ->
        def reads_list = reads instanceof List
            ? is_single_end
                ? reads.findAll { item -> item.name =~ /_singletons.fastq$/ }
                : reads.findAll { item -> item.name =~ /_(1|2)\.fastq$/ }
            : [reads]
        tuple(sample_id, reads_list, is_single_end)
    }
    // repair_read_output_ch.view { it -> "Repair reads output: ${it}" }

    repair_read_size_ch = addSizeTracking(repair_read_size_ch, repair_read_output_ch, "Repair Read Output")

    // Calculate sequence length threshold for repaired reads
    input_to_calculate_sequence_length_threshold_repair = repair_read_output_ch
        .map { sample_id, reads, _is_single_end ->
            tuple(sample_id, reads[0])
        }
        .combine(
            calculate_sequence_length_threshold_script_ch
        )
    calculate_sequence_length_threshold(input_to_calculate_sequence_length_threshold_repair)

    // Quality control process for repaired reads
    repair_read_output_ch
        .join(calculate_sequence_length_threshold.out.length_cutoffs)
        .set { repaired_quality_control_input_ch }

    quality_control(repaired_quality_control_input_ch, 'quality_control', parent_output_dir)
    quality_control_repair_out_ch = quality_control.out.quality_control_output.map { sample_id, reads, is_single_end ->
        def reads_list = reads instanceof List
            ? reads.findAll { item ->
                item.name =~ /_(val_1|val_2)\.fq$/
            }
            : [reads]
        tuple(sample_id, reads_list, is_single_end)
    }
    // Merging quality control outputs and cutoff outputs from both initial(success) and repaired reads
    quality_control_out_ch = quality_control_repair_out_ch
    // .mix(quality_control_repair_out_ch)

    fastqc_process_after_qc(quality_control_repair_out_ch, parent_output_dir, "fastqc_after_quality_control")

    emit:
    quality_control_out_ch
    repair_read_output_ch
    repair_read_size_ch
}
