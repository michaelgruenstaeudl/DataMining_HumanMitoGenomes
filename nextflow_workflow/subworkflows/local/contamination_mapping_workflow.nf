#!/usr/bin/env nextflow

include { contamination_removal } from '../../modules/local/contamination_control.nf'
include { mapping_workflow } from './mapping_workflow.nf'

// Input: tuple(meta, reads), contamination database, confidence threshold,
// and the shared reference channel. Outputs retain tuple(meta, reads).
workflow contamination_mapping_workflow {
    take:
    quality_control_reads_ch
    contamination_db_ch
    contamination_confidence_threshold
    reference_ch

    main:
    contamination_input_ch = quality_control_reads_ch.combine(contamination_db_ch)
    contamination_removal(contamination_input_ch, contamination_confidence_threshold)

    contamination_reads_ch = contamination_removal.out.contamination_removal_fastq_output.map { meta, reads ->
        tuple(meta, reads instanceof List ? reads : [reads])
    }

    mapping_workflow(contamination_reads_ch, reference_ch)
    mapping_reads_ch = mapping_workflow.out.mapping_process_output.map { meta, reads ->
        tuple(meta, reads instanceof List ? reads : [reads])
    }

    emit:
    contamination_reads_ch
    contamination_removal_publish_dir = contamination_removal.out.contamination_removal_publish_dir
    mapping_reads_ch
    index_dir = mapping_workflow.out.index_dir
    mapping_process_publish_dir = mapping_workflow.out.mapping_process_publish_dir
}
