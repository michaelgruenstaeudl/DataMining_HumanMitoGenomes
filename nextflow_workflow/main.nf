#!/usr/bin/env nextflow
include { qualityControl } from './modules/qualityControl.nf'
include { calculate_sequence_length_threshold } from './modules/compute_sequence_length_threshold.nf'
include { mapping_process } from './modules/mapping.nf'
include { novoplast_process } from './modules/novoplasty_assembler.nf'

params.input1 = '/home/b_thapamagar/BioInformatics/NCBIrecordMining/SRA_data/SRR6245224_1.fastq'
params.input2 = '/home/b_thapamagar/BioInformatics/NCBIrecordMining/SRA_data/SRR6245224_2.fastq'
params.reference_fasta = '/home/b_thapamagar/BioInformatics/NCBIrecordMining/nextflow_workflow/reference.fasta'
params.seed_mito = '/home/b_thapamagar/BioInformatics/NCBIrecordMining/nextflow_workflow/Seed_mito.fasta'
params.config_file = '/home/b_thapamagar/BioInformatics/NCBIrecordMining/nextflow_workflow/config.txt'
workflow {
    input_ch = Channel.fromFilePairs("/home/b_thapamagar/BioInformatics/NCBIrecordMining/SRA_data/*_{1,2}.fastq")
        .map { sample_id, read -> tuple(sample_id, read[0], read[1]) }
    // .view()
    reference_ch = channel.fromPath(params.reference_fasta)
    // .view()

    input_to_calculate_sequence_lengt_threshold = input_ch.map { [it[0], it[1]] }
    // .view { "Input to calculate sequence length threshold: ${it}" }
    calculate_sequence_length_threshold(input_to_calculate_sequence_lengt_threshold)

    //Quality control process
    qualityControl_input_ch = input_ch
        .join(calculate_sequence_length_threshold.out.length_cutoffs)
        .view { "Input to quality control: ${it}" }
    qualityControl(qualityControl_input_ch)

    // Mapping process
    mapping_input_ch = qualityControl.out.quality_control_output
        .map { qualityControl_output_ch ->
            def (sample_id, _fastqc1_html, _fastqc1_zip, trimmed_fastq1, _trimmed_fastq1_report, _fastqc2_html, _fastqc2_zip, trimmed_fastq2, _trimmed_fastq2_report) = qualityControl_output_ch
            tuple(sample_id, trimmed_fastq1, trimmed_fastq2)
        }
        .combine(reference_ch)
    // .view { "Input to mapping: ${it}" }

    // combined_ch = sample_id_ch.combine(trimmed_fastq1_ch).combine(trimmed_fastq2_ch)

    // mapping_process(mapping_input_ch)

    // De novo assembly process
    seed_mito_ch = channel.fromPath(params.seed_mito)
    // .view()
    config_file_ch = channel.fromPath(params.config_file)
    // .view()


    denovo_assmebly_input_ch = qualityControl.out.quality_control_output
        .map { qualityControl_output_ch ->
            def (sample_id, _fastqc1_html, _fastqc1_zip, trimmed_fastq1, _trimmed_fastq1_report, _fastqc2_html, _fastqc2_zip, trimmed_fastq2, _trimmed_fastq2_report) = qualityControl_output_ch
            tuple(sample_id, trimmed_fastq1, trimmed_fastq2)
        }
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

    // denovo_assmebly_input_ch = mapping_process.out.mapping_process_output
    //     .map { mapping_process_output_ch ->
    //         def (sample_id, filtered_fq1, filtered_fq2, _sai1, _sai2, _sam) = mapping_process_output_ch

    //         tuple(sample_id, filtered_fq1, filtered_fq2)
    //     }
    //     .join(
    //         calculate_sequence_length_threshold.out.length_cutoffs.map { length_cutoffs_ch ->
    //             def (sample_id, _lower_cutoff, upper_cutoff) = length_cutoffs_ch
    //             def read_length = upper_cutoff + 1
    //             def insert_size = upper_cutoff * 2
    //             tuple(sample_id, read_length, insert_size)
    //         }
    //     )
    //     .combine(seed_mito_ch)
    //     .combine(config_file_ch)
    //     .view { "Input to denovo assembly: ${it}" }


    novoplast_process(denovo_assmebly_input_ch)

    println("Finished all processes")
}
