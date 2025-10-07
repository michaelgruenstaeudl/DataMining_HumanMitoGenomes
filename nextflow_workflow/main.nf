#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { calculate_sequence_length_threshold } from './modules/compute_sequence_length_threshold.nf'
include { mapping_process } from './modules/mapping.nf'
include { novoplast_process } from './modules/novoplasty_assembler.nf'
include { contamination_removal } from "./modules/contamination_control.nf"
include { trackFileSizes ; WriteCSV } from './lib/utils.nf'
include { quality_control_workflow } from './sub_workflows/quality_control_workflow.nf'

workflow {
    input_ch = Channel.fromFilePairs(params.input_files_directory + '/*_{1,2}.fastq')
        .map { sample_id, read -> tuple(sample_id, read[0], read[1]) }
    // .view()

    size_ch = Channel.empty()

    size_ch = size_ch.concat(
        input_ch.map { sample_id, r1, _r2 ->
            trackFileSizes(tuple(sample_id, r1), "Initial Input Files", "Read_1")
        }
    )
    size_ch = size_ch.concat(
        input_ch.map { sample_id, _r1, r2 ->
            trackFileSizes(tuple(sample_id, r2), "Initial Input Files", "Read_2")
        }
    )

    reference_ch = channel.fromPath(params.reference_fasta)
    // .view()
    seed_mito_ch = channel.fromPath(params.seed_mito)
    // .view()
    config_file_ch = channel.fromPath(params.config_file)
    // .view()

    // Calculate sequence length threshold process
    input_to_calculate_sequence_length_threshold = input_ch.map { [it[0], it[1]] }
    // .view { "Input to calculate sequence length threshold: ${it}" }
    calculate_sequence_length_threshold(input_to_calculate_sequence_length_threshold)
    // calculate_sequence_length_threshold.out.length_cutoffs.view { "Output from calculate sequence length threshold: ${it}" }

    //Quality control process
    qualityControl_input_ch = input_ch.join(calculate_sequence_length_threshold.out.length_cutoffs)

    quality_control_workflow(qualityControl_input_ch)
    // quality_control_workflow.out.quality_control_out_ch.view { "Quality Control Output: ${it}" }

    size_ch = size_ch.mix(quality_control_workflow.out.repair_read_size_ch)
    size_ch = size_ch.concat(
        quality_control_workflow.out.quality_control_out_ch.map { sample_id, r1, _r2 ->
            trackFileSizes(tuple(sample_id, r1), "Quality Control Output", "Read_1")
        }
    )
    size_ch = size_ch.concat(
        quality_control_workflow.out.quality_control_out_ch.map { sample_id, _r1, r2 ->
            trackFileSizes(tuple(sample_id, r2), "Quality Control Output", "Read_2")
        }
    )

    //contamination control process
    contamination_db_channel = channel.fromPath(params.contamination_db)
    contamination_control_input_ch = quality_control_workflow.out.quality_control_out_ch.combine(contamination_db_channel)
    contamination_removal(contamination_control_input_ch)

    size_ch = size_ch.concat(
        contamination_removal.out.contamination_removal_fastq_output.map { sample_id, r1, _r2 ->
            trackFileSizes(tuple(sample_id, r1), "Contamination Output", "Read_1")
        }
    )
    size_ch = size_ch.concat(
        contamination_removal.out.contamination_removal_fastq_output.map { sample_id, _r1, r2 ->
            trackFileSizes(tuple(sample_id, r2), "Contamination Output", "Read_2")
        }
    )

    // Mapping process
    mapping_input_ch = contamination_removal.out.contamination_removal_fastq_output.combine(reference_ch)
    // .view { "Input to mapping: ${it}" }
    mapping_process(mapping_input_ch)

    size_ch = size_ch.concat(
        mapping_process.out.mapping_process_output.map { sample_id, r1, _r2 ->
            trackFileSizes(tuple(sample_id, r1), "Mapping Output", "Read_1")
        }
    )
    size_ch = size_ch.concat(
        mapping_process.out.mapping_process_output.map { sample_id, _r1, r2 ->
            trackFileSizes(tuple(sample_id, r2), "Mapping Output", "Read_2")
        }
    )
    csv_lines_ch = size_ch.map { tuple ->
        tuple.join(",")
    }
    // .view { "file sizes: ${it}" }
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
        .view { "Input to denovo assembly: ${it}" }


    novoplast_process(denovo_assmebly_input_ch)

    println("Finished all processes")
}
