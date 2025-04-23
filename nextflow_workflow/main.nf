#!/usr/bin/env nextflow
include { qualityControl } from './modules/qualityControl.nf'
include { calculate_sequence_length_threshold } from './modules/compute_sequence_length_threshold.nf'
params.input1 = '/home/b_thapamagar/BioInformatics/NCBIrecordMining/SRA_data/SRR6245232_1.fastq'
params.input2 = '/home/b_thapamagar/BioInformatics/NCBIrecordMining/SRA_data/SRR6245232_2.fastq'
// 
// println("Input channel: ${input_ch}")
// qualityControl(input_ch)

// println("Finished running the workflow")

workflow {
    // input_ch2 = Channel.fromPath(params.input2, checkIfExists: true).view()

    // input_ch1 = Channel.fromPath(params.input1, checkIfExists: true)
    input_ch = Channel.fromFilePairs("/home/b_thapamagar/BioInformatics/NCBIrecordMining/SRA_data/*_{1,2}.fastq")
        .map { sample_id, read -> tuple(sample_id, read[0], read[1]) }
        .view()
    input_to_calculate_sequence_lengt_threshold = input_ch.map { [it[0], it[1]] }.view { "Input to calculate sequence length threshold: ${it}" }
    calculate_sequence_length_threshold(input_to_calculate_sequence_lengt_threshold)
    println(calculate_sequence_length_threshold.out.length_cutoffs.map { it[0] }.view())
    qualityControl(input_ch, calculate_sequence_length_threshold.out.length_cutoffs)
}
