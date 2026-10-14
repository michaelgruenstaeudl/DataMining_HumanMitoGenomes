#!/usr/bin/env nextflow

include { generate_sam } from '../../modules/local/generate_sam.nf'
include { sam_to_bam } from '../../modules/local/sam_to_bam.nf'
include { calculate_evenness } from '../../modules/local/calculate_evenness.nf'
include { repair_reads } from '../../modules/local/repair_reads.nf'

workflow evenness_calculation_workflow {
    take:
    evenness_calculation_input_channel
    phase

    main:
    fasta_file_ch = evenness_calculation_input_channel.map { meta, _fastq_list, _gb_file, fasta_file ->
        tuple(meta, fasta_file)
    }


    gb_file_ch = evenness_calculation_input_channel.map { meta, _fastq_list, gb_file, _fasta_file ->
        tuple(meta, gb_file)
    }
    generate_sam_input_ch = evenness_calculation_input_channel
        .map { meta, fastq_list, _gb_file, fasta_file ->
            tuple(meta, fastq_list, fasta_file)
        }
    // }

    generate_sam(
        generate_sam_input_ch,
        phase,
    )
    sam_to_bam_input_ch = fasta_file_ch.join(generate_sam.out.sam_output_ch)
    sam_to_bam(
        sam_to_bam_input_ch,
        phase,
    )
    calculate_evenness_input_ch = gb_file_ch.join(sam_to_bam.out.bam_output_ch)
    calculate_evenness(
        calculate_evenness_input_ch,
        phase,
    )
    evenness_output_ch = calculate_evenness.out.evenness_output_ch

    emit:
    sam_output_publish_dir = generate_sam.out.sam_output_publish_dir
    bam_output_publish_dir = sam_to_bam.out.bam_output_publish_dir
    evenness_output_publish_dir = calculate_evenness.out.evenness_output_publish_dir

    evenness_output_ch
}
