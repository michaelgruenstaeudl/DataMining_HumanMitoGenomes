#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { calculate_sequence_length_threshold } from './modules/compute_sequence_length_threshold.nf'
include { mapping_process } from './modules/mapping.nf'
include { novoplast_process } from './modules/novoplasty_assembler.nf'
include { contamination_removal } from "./modules/contamination_control.nf"
include { addSizeTracking ; WriteCSV } from './lib/utils.nf'
include { quality_control_workflow } from './sub_workflows/quality_control_workflow.nf'

workflow {
    input_ch = Channel.fromFilePairs(params.input_files_directory + '/*_{1,2}.fastq')
        .map { sample_id, read -> tuple(sample_id, read[0], read[1]) }
    // .view()

    size_ch = Channel.empty()

    size_ch = addSizeTracking(size_ch, input_ch, "Initial Input Files")

    reference_ch = channel.fromPath(params.reference_fasta)
    seed_mito_ch = channel.fromPath(params.seed_mito)
    config_file_ch = channel.fromPath(params.config_file)

    // Calculate sequence length threshold process
    input_to_calculate_sequence_length_threshold = input_ch.map { [it[0], it[1]] }
    calculate_sequence_length_threshold(input_to_calculate_sequence_length_threshold)


    //Quality control process
    qualityControl_input_ch = input_ch.join(calculate_sequence_length_threshold.out.length_cutoffs)

    quality_control_workflow(qualityControl_input_ch)
    size_ch = size_ch.mix(quality_control_workflow.out.repair_read_size_ch)
    size_ch = addSizeTracking(size_ch, quality_control_workflow.out.quality_control_out_ch, "Quality Control Output")

    //contamination control process
    contamination_db_channel = channel.fromPath(params.contamination_db)
    contamination_control_input_ch = quality_control_workflow.out.quality_control_out_ch.combine(contamination_db_channel)
    contamination_removal(contamination_control_input_ch)
    size_ch = addSizeTracking(size_ch, contamination_removal.out.contamination_removal_fastq_output, "Contamination Output")


    // Mapping process
    mapping_input_ch = contamination_removal.out.contamination_removal_fastq_output.combine(reference_ch)
    mapping_process(mapping_input_ch)
    size_ch = addSizeTracking(size_ch, mapping_process.out.mapping_process_output, "Mapping Output")

    csv_lines_ch = size_ch.map { tuple ->
        tuple.join(",")
    }
    tmp_csv = csv_lines_ch.collectFile(name: "file_sizes.tmp", newLine: true)

    WriteCSV(tmp_csv)
    size_ch.count().view { "Total tuples in channel: ${it}" }

    // De novo assembly process
    denovo_assmebly_input_ch = mapping_process.out.mapping_process_output
        .join(
            calculate_sequence_length_threshold.out.length_cutoffs.map { length_cutoffs_ch ->
                def (sample_id, _lower_cutoff, upper_cutoff) = length_cutoffs_ch
                def read_length = upper_cutoff.toInteger() + 1
                def insert_size = upper_cutoff.toInteger() * 2
                tuple(sample_id, read_length, insert_size)
            }
        )
        .combine(seed_mito_ch)
        .combine(config_file_ch)


    novoplast_process(denovo_assmebly_input_ch)

    workflow.onComplete {
        println("✅ Finished all processes!")
    }
}
