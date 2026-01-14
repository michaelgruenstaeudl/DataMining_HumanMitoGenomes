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

workflow {
    parent_output_dir = params.outdir

    // Call the read_csv_workflow with params
    read_csv_workflow(params)
    evenness_calculation_workflow_initial_input_ch = read_csv_workflow.out.evenness_calculation_workflow_initial_input_ch
    accession_ch = read_csv_workflow.out.accession_ch
    input_ch = read_csv_workflow.out.input_ch
    gb_channel = read_csv_workflow.out.gb_channel
    fasta_channel = read_csv_workflow.out.fasta_channel

    merge_pacvr_evenness_fig_script_ch = channel.fromPath(params.merge_pacvr_evenness_figure_script)

    evenness_calculation_workflow_initial(evenness_calculation_workflow_initial_input_ch, 'initial_phase', parent_output_dir)
    initial_evenness_output_ch = evenness_calculation_workflow_initial.out.evenness_output_ch

    size_ch = channel.empty()

    size_ch = addSizeTracking(size_ch, input_ch, "Initial Input Files")

    calculate_sequence_length_threshold_script_ch = channel.fromPath(params.calculate_sequence_length_threshold_script)

    //Quality control process
    quality_control_workflow(input_ch, calculate_sequence_length_threshold_script_ch, parent_output_dir)
    size_ch = size_ch.mix(quality_control_workflow.out.repair_read_size_ch)
    size_ch = addSizeTracking(size_ch, quality_control_workflow.out.quality_control_out_ch, "Quality Control Output")

    //Evenness calculation workflow
    final_evenness_calculation_input_ch = quality_control_workflow.out.quality_control_out_ch.join(gb_channel).join(fasta_channel)
    evenness_calculation_workflow_final(final_evenness_calculation_input_ch, "mapping_phase", parent_output_dir)
    contamination_removed_evenness_output_ch = evenness_calculation_workflow_final.out.evenness_output_ch

    merge_evenness_figure_input_ch = accession_ch
        .join(initial_evenness_output_ch)
        .join(contamination_removed_evenness_output_ch)
        .combine(merge_pacvr_evenness_fig_script_ch)

    merge_pacvr_evenness_figures(
        merge_evenness_figure_input_ch,
        parent_output_dir,
    )

    // Write file sizes to CSV
    file_size_csv_lines_ch = size_ch.map { tuple ->
        tuple.join(",")
    }
    file_size_tmp_csv = file_size_csv_lines_ch.collectFile(name: "file_sizes.tmp", newLine: true)

    write_file_sizes_csv(file_size_tmp_csv, true, "file_sizes.csv", parent_output_dir)
    size_ch.count().view { it -> "Total tuples in channel: ${it}" }


    plot_size_diff_input_ch = write_file_sizes_csv.out.csv_output
        .combine(channel.fromPath(parent_output_dir))
        .combine(channel.fromPath(params.visualize_file_size_changes_script))
    plot_size_diff(plot_size_diff_input_ch)


    // MitoZ assembly process
    denovo_assmebly_input_ch = quality_control_workflow.out.quality_control_out_ch
    mitoz_assembler(denovo_assmebly_input_ch, parent_output_dir)
    extract_assembled_genome_metadata_info(mitoz_assembler.out.mitoz_assembler_output, parent_output_dir)

    // Write assembled genome metadata info into CSV
    assembled_genome_csv_lines_ch = extract_assembled_genome_metadata_info.out.metadata_info_out_ch.map { tuple ->
        tuple.join(",")
    }
    tmp_csv = assembled_genome_csv_lines_ch.collectFile(name: "assembled_genome_metadata.tmp", newLine: true)
    write_assembled_metadata_csv(tmp_csv, false, "assembled_genome_metadata.csv", parent_output_dir)


    workflow.onComplete {
        println('✅ Finished all processes!')
    }
}
