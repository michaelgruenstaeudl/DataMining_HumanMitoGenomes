#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { calculate_sequence_length_threshold } from './modules/compute_sequence_length_threshold.nf'
include { mapping_process } from './modules/mapping_process.nf'
include { novoplast_process } from './modules/novoplasty_process.nf'
include { contamination_removal } from './modules/contamination_control.nf'
include { addSizeTracking ; WriteCSV } from './lib/utils.nf'
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
workflow {
    parent_output_dir = params.outdir

    // Call the read_csv_workflow with params
    read_csv_workflow(params)
    // evenness_calculation_workflow_initial_input_ch = read_csv_workflow.out.evenness_calculation_workflow_initial_input_ch
    // accession_ch = read_csv_workflow.out.accession_ch
    input_ch = read_csv_workflow.out.input_ch
    // gb_channel = read_csv_workflow.out.gb_channel
    // fasta_channel = read_csv_workflow.out.fasta_channel

    // merge_pacvr_evenness_fig_script_ch = channel.fromPath(params.merge_pacvr_evenness_figure_script)

    // evenness_calculation_workflow_initial(evenness_calculation_workflow_initial_input_ch, 'initial_phase', parent_output_dir)
    // initial_evenness_output_ch = evenness_calculation_workflow_initial.out.evenness_output_ch

    size_ch = channel.empty()

    size_ch = addSizeTracking(size_ch, input_ch, "Initial Input Files")
    mitoz_assembler(input_ch, parent_output_dir)

    workflow.onComplete {
        println('✅ Finished all processes!')
    }
}
