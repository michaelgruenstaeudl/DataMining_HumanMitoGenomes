#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { calculate_sequence_length_threshold } from './modules/compute_sequence_length_threshold.nf'
include { mapping_process } from './modules/mapping_process.nf'
include { novoplast_process } from './modules/novoplasty_process.nf'
include { contamination_removal } from './modules/contamination_control.nf'
include { addSizeTracking } from './lib/utils.nf'
include {
    write_csv as write_file_sizes_csv ;
    write_csv as write_assembled_metadata_csv ;
    write_csv as write_rotate_genome_csv ;
    write_csv as write_rotate_genome_csv_for_reporting
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
include { rotate_circular_genome_workflow } from './sub_workflows/rotate_circular_genome_workflow.nf'
include { normalize_complete_genome_length } from './sub_workflows/normalize_complete_genome_length.nf'
workflow {
    parent_output_dir = params.outdir
    contamination_confidence_threshold = params.contamination_confidence_threshold
    config_file_ch = channel.fromPath(params.config_file)
    // Call the read_csv_workflow with params
    read_csv_workflow(params)
    // evenness_calculation_workflow_initial_input_ch = read_csv_workflow.out.evenness_calculation_workflow_initial_input_ch
    accession_ch = read_csv_workflow.out.accession_ch
    input_ch = read_csv_workflow.out.input_ch
    gb_channel = read_csv_workflow.out.gb_channel
    fasta_channel = read_csv_workflow.out.fasta_channel

    size_ch = channel.empty()

    size_ch = addSizeTracking(size_ch, input_ch, "Initial Input Files")
    // Created once here and passed to mapping workflow, which will then pass it to the mapping processes. 
    //This way we avoid creating multiple index files if the same reference is used for multiple samples.
    reference_ch = channel.value(file(params.reference_fasta))

    calculate_sequence_length_threshold_script_ch = channel.fromPath(params.calculate_sequence_length_threshold_script)
    rotation_seed_fasta_ch = channel.fromPath(params.rotation_seed_fasta)
    trim_mitogenome_duplicate_script_ch = channel.fromPath(params.trim_mitogenome_duplicate_script)
    //Quality control process
    quality_control_workflow(
        input_ch,
        calculate_sequence_length_threshold_script_ch,
        parent_output_dir,
    )
    cutoffs_ch = quality_control_workflow.out.cutoffs_ch
    evenness_calculation_workflow_initial_input_ch = quality_control_workflow.out.repair_read_output_ch
        .join(gb_channel)
        .join(fasta_channel)
    // quality_control_workflow.out.quality_control_out_ch.view { it -> "QC output: ${it}" }
    size_ch = size_ch.mix(quality_control_workflow.out.repair_read_size_ch)
    size_ch = addSizeTracking(size_ch, quality_control_workflow.out.quality_control_out_ch, "Quality Control Output")

    //contamination control process
    contamination_db_channel = channel.fromPath(params.contamination_db)
    contamination_control_input_ch = quality_control_workflow.out.quality_control_out_ch.combine(contamination_db_channel)
    contamination_removal(
        contamination_control_input_ch,
        contamination_confidence_threshold,
        parent_output_dir,
    )
    // contamination_removal_out_ch = contamination_removal.out.contamination_removal_fastq_output
    contamination_removal_out_ch = contamination_removal.out.contamination_removal_fastq_output.map { sample_id, reads, is_single_end ->
        def reads_list = is_single_end
            ? [reads]
            : reads
        tuple(sample_id, reads_list, is_single_end)
    }
    // contamination_removal_out_ch.view { it -> "Contamination removal output: ${it}" }
    size_ch = addSizeTracking(size_ch, contamination_removal_out_ch, "Contamination Output")

    // Mapping process
    // index_files_ch = channel.fromPath("${params.index_directory}/*").collect()
    // mapping_input_ch = contamination_removal_out_ch
    //.combine(index_files_ch)
    mapping_workflow(contamination_removal_out_ch, reference_ch, parent_output_dir)
    mapping_process_output = mapping_workflow.out.mapping_process_output.map { sample_id, reads, is_single_end ->
        def reads_list = is_single_end
            ? [reads]
            : reads
        tuple(sample_id, reads_list, is_single_end)
    }

    size_ch = addSizeTracking(size_ch, mapping_process_output, "Mapping Output")

    if (params.is_calculate_evenness) {
        //Evenness calculation workflow
        merge_pacvr_evenness_fig_script_ch = channel.fromPath(params.merge_pacvr_evenness_figure_script)

        //evenness calculation for original reads before cleaning
        evenness_calculation_workflow_initial(
            evenness_calculation_workflow_initial_input_ch,
            'initial_phase',
            parent_output_dir,
        )
        initial_evenness_output_ch = evenness_calculation_workflow_initial.out.evenness_output_ch

        // evenness calculation for cleaned reads after mapping
        final_evenness_calculation_input_ch = mapping_process_output
            .join(gb_channel)
            .join(fasta_channel)
        evenness_calculation_workflow_final(
            final_evenness_calculation_input_ch,
            "mapping_phase",
            parent_output_dir,
        )
        cleaned_reads_evenness_output_ch = evenness_calculation_workflow_final.out.evenness_output_ch

        merge_evenness_figure_input_ch = accession_ch
            .join(initial_evenness_output_ch)
            .join(cleaned_reads_evenness_output_ch)
            .combine(merge_pacvr_evenness_fig_script_ch)

        merge_pacvr_evenness_figures(
            merge_evenness_figure_input_ch,
            parent_output_dir,
        )
    }

    // Write file sizes to CSV
    file_size_csv_lines_ch = size_ch.map { tuple ->
        tuple.join(",")
    }
    file_size_tmp_csv = file_size_csv_lines_ch.collectFile(name: "file_sizes.tmp", newLine: true)

    write_file_sizes_csv(file_size_tmp_csv, "file_sizes", "file_sizes.csv", parent_output_dir)
    size_ch.count().view { it -> "Total tuples in channel: ${it}" }


    plot_size_diff_input_ch = write_file_sizes_csv.out.csv_output
        .combine(channel.fromPath(parent_output_dir))
        .combine(channel.fromPath(params.visualize_file_size_changes_script))
    plot_size_diff(plot_size_diff_input_ch)


    // MitoZ assembly process
    mitoz_assembler(mapping_process_output, parent_output_dir)

    extract_assembled_genome_metadata_info(
        mitoz_assembler.out.mitoz_assembler_output,
        parent_output_dir,
    )

    // Write assembled genome metadata info into CSV
    assembled_genome_csv_lines_ch = extract_assembled_genome_metadata_info.out.metadata_info_out_ch.map { tuple ->
        tuple.join(",")
    }
    tmp_csv = assembled_genome_csv_lines_ch.collectFile(
        name: "assembled_genome_metadata.tmp",
        newLine: true,
    )
    write_assembled_metadata_csv(
        tmp_csv,
        "assembled_genome_metadata",
        "assembled_genome_metadata.csv",
        parent_output_dir,
    )

    normalize_fasta_input_ch = mapping_process_output
        .join(
            cutoffs_ch.map { length_cutoffs_ch ->
                def (sample_id, _lower_cutoff, upper_cutoff) = length_cutoffs_ch
                def read_length = upper_cutoff.toInteger() + 1
                def insert_size = upper_cutoff.toInteger() * 2
                tuple(sample_id, read_length, insert_size)
            }
        )
        .join(mitoz_assembler.out.mitoz_assembler_output)
        .combine(config_file_ch)
        .combine(trim_mitogenome_duplicate_script_ch)
    normalize_complete_genome_length(
        normalize_fasta_input_ch,
        trim_mitogenome_duplicate_script_ch,
        parent_output_dir,
    )

    //Circular genome rotation process
    rotate_genome_input_ch = normalize_complete_genome_length.out.normalized_fasta_ch
        .join(fasta_channel)
        .combine(rotation_seed_fasta_ch)
        .combine(channel.value(params.rotation_mismatch_threshold))

    rotate_circular_genome_workflow(rotate_genome_input_ch, parent_output_dir)

    rotated_assembled_ch = rotate_circular_genome_workflow.out.rotated_assembled_ch
    rotated_official_ch = rotate_circular_genome_workflow.out.rotated_official_ch


    rotated_fasta_ch = rotated_assembled_ch
        .join(rotated_official_ch)
        .map { tuple ->
            tuple.join(",")
        }
    rotated_fasta_tmp_csv = rotated_fasta_ch.collectFile(name: "rotated_genome_info.tmp", newLine: true)

    write_rotate_genome_csv(
        rotated_fasta_tmp_csv,
        "rotated_genome_info",
        "rotated_genome_info.csv",
        parent_output_dir,
    )

    // For reporting, we want to include the paths to the rotated genome files instead of the actual sequences. 
    // So we construct new lines for the CSV that have the sample ID and the paths to the rotated assembled 
    // and official genome files, and then write that to a separate CSV for reporting purposes.
    rotated_fasta_ch_for_reporting = rotated_assembled_ch
        .join(rotated_official_ch)
        .map { sample_id, _assembled, _official ->
            // Construct results dir paths — NOT work dir paths
            def assembled_result = "${params.outdir}/${sample_id}/rotated_output/assembled_genome/rotated_${sample_id}.fasta"
            def official_result = "${params.outdir}/${sample_id}/rotated_output/official_genome/rotated_${sample_id}.fasta"

            "${sample_id},${assembled_result},${official_result}"
        }

    rotated_fasta_reporting_tmp_csv = rotated_fasta_ch_for_reporting.collectFile(name: "rotated_genome_reporting_info.tmp", newLine: true)

    write_rotate_genome_csv_for_reporting(
        rotated_fasta_reporting_tmp_csv,
        "rotated_genome_info",
        "rotated_genome_reporting_info.csv",
        parent_output_dir,
    )

    workflow.onComplete {
        println('✅ Finished all processes!')
    }
}
