#!/usr/bin/env nextflow

include {
    write_csv as write_cutoff_info_csv;
    write_csv as write_file_sizes_csv;
    write_csv as write_assembled_metadata_csv;
    write_csv as write_normalized_genome_metadata_csv;
    write_csv as write_rotate_genome_csv;
    write_csv as write_rotate_genome_csv_for_reporting;
    write_csv as write_assembled_haplotype_info_csv;
    write_csv as write_official_haplotype_info_csv
} from '../../modules/local/write_csv.nf'
include { merge_pacvr_evenness_figures } from '../../modules/local/merge_pacvr_evenness_figures.nf'
include {
    plot_size_diff;
    plot_size_diff as plot_genome_comparison
} from '../../modules/local/plot_size_diff.nf'

// Report generation consumes completed sample-level results and emits only
// report files and report-related publish channels used by the top-level workflow.
workflow report_generation_workflow {
    take:
    cutoffs_ch
    size_ch
    assembled_metadata_info_ch
    normalized_metadata_info_ch
    rotated_genome_ch
    rotated_assembled_ch
    rotated_official_ch
    assembled_haplotype_info_ch
    official_haplotype_info_ch
    calculate_evenness
    initial_evenness_output_ch
    final_evenness_output_ch
    merge_pacvr_evenness_fig_script_ch
    visualize_file_size_changes_script_ch
    genome_comparison_plot_script_ch

    main:
    cutoff_lines_ch = cutoffs_ch.map { meta, lower_cutoff, upper_cutoff ->
        [meta.id, lower_cutoff, upper_cutoff].join(',')
    }
    cutoff_tmp_csv = cutoff_lines_ch.collectFile(name: 'cutoff_info.tmp', newLine: true)
    write_cutoff_info_csv(cutoff_tmp_csv, 'cutoff_info', 'cutoff_info.csv')

    file_size_csv_lines_ch = size_ch.map { tuple -> tuple.join(',') }
    file_size_tmp_csv = file_size_csv_lines_ch.collectFile(name: 'file_sizes.tmp', newLine: true)
    write_file_sizes_csv(file_size_tmp_csv, 'file_sizes', 'file_sizes.csv')

    plot_size_diff_input_ch = write_file_sizes_csv.out.csv_output
        .combine(channel.fromPath(workflow.outputDir))
        .combine(visualize_file_size_changes_script_ch)
    plot_size_diff(plot_size_diff_input_ch)

    assembled_genome_csv_lines_ch = assembled_metadata_info_ch.map { meta, topology, sequence_length ->
        [meta.id, topology, sequence_length].join(',')
    }
    assembled_genome_csv_tmp = assembled_genome_csv_lines_ch.collectFile(name: 'assembled_genome_metadata.tmp', newLine: true)
    write_assembled_metadata_csv(assembled_genome_csv_tmp, 'assembled_genome_metadata', 'assembled_genome_metadata.csv')

    normalized_genome_csv_lines_ch = normalized_metadata_info_ch.map { meta, topology, sequence_length ->
        [meta.id, topology, sequence_length].join(',')
    }
    normalized_genome_csv_tmp = normalized_genome_csv_lines_ch.collectFile(name: 'normalized_genome_metadata.tmp', newLine: true)
    write_normalized_genome_metadata_csv(normalized_genome_csv_tmp, 'normalized_genome_metadata', 'normalized_genome_metadata.csv')

    rotated_fasta_ch = rotated_genome_ch.map { meta, assembled_fasta, official_fasta ->
        [meta.id, assembled_fasta, official_fasta].join(',')
    }
    rotated_fasta_tmp_csv = rotated_fasta_ch.collectFile(name: 'rotated_genome_info.tmp', newLine: true)
    write_rotate_genome_csv(rotated_fasta_tmp_csv, 'rotated_genome_info', 'rotated_genome_info.csv')

    rotated_fasta_ch_for_reporting = rotated_assembled_ch
        .join(rotated_official_ch)
        .map { meta, _assembled, _official ->
            def assembled_result = "${workflow.outputDir}/${meta.id}/rotated_output/assembled_genome/rotated_${meta.id}.fasta"
            def official_result = "${workflow.outputDir}/${meta.id}/rotated_output/official_genome/rotated_${meta.id}.fasta"
            "${meta.id},${assembled_result},${official_result}"
        }
    rotated_fasta_reporting_tmp_csv = rotated_fasta_ch_for_reporting.collectFile(name: 'rotated_genome_reporting_info.tmp', newLine: true)
    write_rotate_genome_csv_for_reporting(rotated_fasta_reporting_tmp_csv, 'rotated_genome_info', 'rotated_genome_reporting_info.csv')

    plot_genome_comparison_input_ch = write_rotate_genome_csv.out.csv_output
        .combine(channel.fromPath(workflow.outputDir))
        .combine(genome_comparison_plot_script_ch)
    plot_genome_comparison(plot_genome_comparison_input_ch)

    assembled_haplotype_csv_lines_ch = assembled_haplotype_info_ch.map { meta, haplotype, quality ->
        [meta.id, haplotype, quality].join(',')
    }
    assembled_haplotype_tmp_csv = assembled_haplotype_csv_lines_ch.collectFile(name: 'assembled_haplotype_info.tmp', newLine: true)
    write_assembled_haplotype_info_csv(assembled_haplotype_tmp_csv, 'haplotype_info', 'assembled_haplotype_info.csv')

    official_haplotype_csv_lines_ch = official_haplotype_info_ch.map { meta, haplotype, quality ->
        [meta.id, haplotype, quality].join(',')
    }
    official_haplotype_tmp_csv = official_haplotype_csv_lines_ch.collectFile(name: 'official_haplotype_info.tmp', newLine: true)
    write_official_haplotype_info_csv(official_haplotype_tmp_csv, 'haplotype_info', 'official_haplotype_info.csv')

    if (calculate_evenness) {
        merge_evenness_figure_input_ch = initial_evenness_output_ch
            .join(final_evenness_output_ch)
            .combine(merge_pacvr_evenness_fig_script_ch)
            .map { meta, initial_figure, final_figure, merge_script ->
                tuple(meta, initial_figure, final_figure, merge_script)
        }
        merge_pacvr_evenness_figures(merge_evenness_figure_input_ch)
        merged_pacvr_figures_output_publish_dir = merge_pacvr_evenness_figures.out.merged_pacvr_figures_output_ch
        merged_pacvr_figures_files = merge_pacvr_evenness_figures.out.merged_pacvr_figures_files
    }
    else {
        merged_pacvr_figures_output_publish_dir = channel.empty()
        merged_pacvr_figures_files = channel.empty()
    }

    emit:
    cutoff_info_csv = write_cutoff_info_csv.out.csv_output
    file_size_csv = write_file_sizes_csv.out.csv_output
    assembled_genome_metadata_csv = write_assembled_metadata_csv.out.csv_output
    normalized_genome_metadata_csv = write_normalized_genome_metadata_csv.out.csv_output
    rotated_genome_metadata_csv = write_rotate_genome_csv.out.csv_output
    rotated_genome_csv_for_reporting = write_rotate_genome_csv_for_reporting.out.csv_output
    assembled_haplotype_info_csv = write_assembled_haplotype_info_csv.out.csv_output
    official_haplotype_info_csv = write_official_haplotype_info_csv.out.csv_output
    merged_pacvr_figures_output_publish_dir
    merged_pacvr_figures_files
}
