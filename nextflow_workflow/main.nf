#!/usr/bin/env nextflow
include { qualityControl } from './modules/qualityControl.nf'
include { calculate_sequence_length_threshold } from './modules/compute_sequence_length_threshold.nf'
include { mapping_process } from './modules/mapping.nf'
params.input1 = '/home/b_thapamagar/BioInformatics/NCBIrecordMining/SRA_data/SRR6245224_1.fastq'
params.input2 = '/home/b_thapamagar/BioInformatics/NCBIrecordMining/SRA_data/SRR6245224_2.fastq'
params.reference_fasta = '/home/b_thapamagar/BioInformatics/NCBIrecordMining/SRA_data/reference.fasta'

workflow {
    input_ch = Channel.fromFilePairs("/home/b_thapamagar/BioInformatics/NCBIrecordMining/SRA_data/*_{1,2}.fastq")
        .map { sample_id, read -> tuple(sample_id, read[0], read[1]) }
        .view()
    reference_ch = channel
        .fromPath(params.reference_fasta)
        .view()
    input_to_calculate_sequence_lengt_threshold = input_ch.map { [it[0], it[1]] }.view { "Input to calculate sequence length threshold: ${it}" }
    calculate_sequence_length_threshold(input_to_calculate_sequence_lengt_threshold)
    println(calculate_sequence_length_threshold.out.length_cutoffs.view())
    qualityControl(input_ch, calculate_sequence_length_threshold.out.length_cutoffs)

    // Collect the outputs from qualityControl
    sample_id_ch = qualityControl.out.sample_id.collect()
    trimmed_fastq1_ch = qualityControl.out.trimmed_fastq1.collect()
    trimmed_fastq2_ch = qualityControl.out.trimmed_fastq2.collect()

    combined_ch = sample_id_ch.combine(trimmed_fastq1_ch).combine(trimmed_fastq2_ch)

    mapping_process(combined_ch, reference_ch)
}
