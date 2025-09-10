#!/usr/bin/env nextflow
include { qualityControl } from './modules/qualityControl.nf'
include { calculate_sequence_length_threshold } from './modules/compute_sequence_length_threshold.nf'
include { mapping_process } from './modules/mapping.nf'
include { novoplast_process } from './modules/novoplasty_assembler.nf'
include { contamination_removal } from "./modules/contamination_control.nf"

params.input1 = '/home/b_thapamagar/BioInformatics/NCBIrecordMining/SRA_data/SRR6245224_1.fastq'
params.input2 = '/home/b_thapamagar/BioInformatics/NCBIrecordMining/SRA_data/SRR6245224_2.fastq'

workflow {
    input_ch = Channel.fromFilePairs(params.input_files_directory + '/*_{1,2}.fastq')
        .map { sample_id, read -> tuple(sample_id, read[0], read[1]) }
    // .view()
    reference_ch = channel.fromPath(params.reference_fasta)
    // .view()


    // Calculate sequence length threshold process
    input_to_calculate_sequence_length_threshold = input_ch.map { [it[0], it[1]] }
    // .view { "Input to calculate sequence length threshold: ${it}" }
    calculate_sequence_length_threshold(input_to_calculate_sequence_length_threshold)
    calculate_sequence_length_threshold.out.length_cutoffs.view { "Output from calculate sequence length threshold: ${it}" }


    //Quality control process
    qualityControl_input_ch = input_ch
        .join(calculate_sequence_length_threshold.out.length_cutoffs)
        .view { "Input to quality control: ${it}" }
    qualityControl(qualityControl_input_ch)


    //contamination control process
    contamination_db_channel = channel.fromPath(params.contamination_db)
    contamination_control_input_ch = qualityControl.out.quality_control_output.combine(contamination_db_channel)
    contamination_removal(contamination_control_input_ch)


    // Mapping process
    mapping_input_ch = contamination_removal.out.contamination_removal_fastq_output.combine(reference_ch)
    // .view { "Input to mapping: ${it}" }
    mapping_process(mapping_input_ch)


    // De novo assembly process
    seed_mito_ch = channel.fromPath(params.seed_mito)
    // .view()
    config_file_ch = channel.fromPath(params.config_file)
    // .view()

    denovo_assmebly_input_ch = mapping_process.out.mapping_process_output
        .join(
            calculate_sequence_length_threshold.out.length_cutoffs.map { length_cutoffs_ch ->
                def (sample_id, _lower_cutoff, upper_cutoff) = length_cutoffs_ch
                def read_length = upper_cutoff + 1
                def insert_size = upper_cutoff * 2
                tuple(sample_id, read_length, insert_size)
            }
        )
        .combine(seed_mito_ch)
        .combine(config_file_ch)
        .view { "Input to denovo assembly: ${it}" }


    novoplast_process(denovo_assmebly_input_ch)

    println("Finished all processes")
}
