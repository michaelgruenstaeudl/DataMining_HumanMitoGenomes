#!/usr/bin/env nextflow

include { mapping_process } from '../modules/mapping_process.nf'
include { bbmap_workflow } from './bbmap_workflow.nf'
include { mapping_bwa } from '../modules/mapping_bwa.nf'
include { extract_mapped_fastq } from '../modules/extract_mapped_fastq.nf'

workflow mapping_workflow {
    take:
    mapping_input_ch
    reference_channel
    parent_output_dir

    main:
    // mapping_input_ch.view { it -> "Mapping input: ${it}" }

    //BWA_MEM2 mapping process
    // mapping_bwa(mapping_input_ch, parent_output_dir)
    // extract_mapped_fastq(mapping_bwa.out.mapping_bwa_output, parent_output_dir)
    // mapping_process_output = extract_mapped_fastq.out.extract_mapped_fastq_output_ch

    // FASTQSIFTER mapping process
    // mapping_process(input_ch, parent_output_dir)
    // mapping_process_output = mapping_process.out.mapping_process_output


    //####BBMAP mapping process####
    bbmap_workflow(mapping_input_ch, reference_channel, parent_output_dir)
    mapping_process_output = bbmap_workflow.out.mapping_process_output

    emit:
    mapping_process_output
}
