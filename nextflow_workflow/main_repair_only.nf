#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { calculate_sequence_length_threshold } from './modules/compute_sequence_length_threshold.nf'
include { mapping_process } from './modules/mapping_process.nf'
include { novoplast_process } from './modules/novoplasty_process.nf'
include { contamination_removal } from './modules/contamination_control.nf'
include { addSizeTracking } from './lib/utils.nf'
include {
    write_csv as write_file_sizes_csv ;
    write_csv as write_assembled_metadata_csv
} from './modules/write_csv.nf'
include { quality_control_workflow } from './sub_workflows/quality_control_workflow.nf'
include { extract_seed } from './modules/extract_seed.nf'
include { merge_pacvr_evenness_figures } from './modules/merge_pacvr_evenness_figures.nf'
include { read_csv_workflow } from './sub_workflows/read_csv_workflow.nf'
include {
    evenness_calculation_workflow as evenness_calculation_workflow_initial ;
    evenness_calculation_workflow as evenness_calculation_workflow_final
} from './sub_workflows/evenness_calculation_workflow.nf'
include { plot_size_diff } from './modules/plot_size_diff.nf'
include { mapping_workflow } from './sub_workflows/mapping_workflow.nf'
include { mapping_bbmap } from './modules/mapping_bbmap.nf'
include { mitoz_assembler } from './modules/mitoz_assembler.nf'
include { extract_assembled_genome_metadata_info } from './modules/extract_assembled_genome_metadata_info.nf'
include { repair_reads } from './modules/repair_reads.nf'

workflow {
    parent_output_dir = params.outdir
    // Call the read_csv_workflow with params
    read_csv_workflow(params)
    input_ch = read_csv_workflow.out.input_ch

    size_ch = channel.empty()

    size_ch = addSizeTracking(size_ch, input_ch, "Initial Input Files")
    // Created once here and passed to mapping workflow, which will then pass it to the mapping processes. 
    //This way we avoid creating multiple index files if the same reference is used for multiple samples.
    repair_reads(input_ch, 'pre_mapping', parent_output_dir)

    workflow.onComplete {
        println('✅ Finished all processes!')
    }
}
