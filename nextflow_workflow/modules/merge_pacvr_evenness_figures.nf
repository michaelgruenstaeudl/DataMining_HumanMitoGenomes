#!/usr/bin/env nextflow

process merge_pacvr_evenness_figures {

    tag "${sample_id}"

    // container "texlive/texlive:latest"

    publishDir "${parent_output_dir}/${sample_id}", mode: 'copy', overwrite: true
    // 2️⃣ Also collect specific outputs in one shared folder
    publishDir "${parent_output_dir}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), val(accession_number), path(pacvr_evenness_initial_figure), path(pacvr_evenness_figure_after_contamination_removal), path(merge_pacvr_evenness_fig_script)
    val parent_output_dir

    output:
    // tuple val(sample_id), path("pacvr/${sample_id}_original_CoverageViz.pdf"), emit: evenness_output_ch
    path "combined_coverage_visualization/*"

    script:

    """
    mkdir -p combined_coverage_visualization
    bash ${merge_pacvr_evenness_fig_script} ${sample_id} ${accession_number} ${pacvr_evenness_initial_figure} ${pacvr_evenness_figure_after_contamination_removal}
    """
}
