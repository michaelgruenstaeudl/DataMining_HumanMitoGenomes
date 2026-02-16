#!/usr/bin/env nextflow

include { mapping_process } from '../modules/mapping_process.nf'
include { mapping_bbmap } from '../modules/mapping_bbmap.nf'
include { mapping_bwa } from '../modules/mapping_bwa.nf'
include { extract_mapped_fastq } from '../modules/extract_mapped_fastq.nf'
include { bbmap_index } from '../modules/bbmap_index.nf'

workflow bbmap_workflow {
    take:
    mapping_input_ch
    reference_channel
    parent_output_dir

    main:

    index_dir = bbmap_index(reference_channel, parent_output_dir)
    mapping_bbmap(mapping_input_ch, index_dir, parent_output_dir)
    mapping_process_output = mapping_bbmap.out.mapping_process_output

    emit:
    mapping_process_output
}
