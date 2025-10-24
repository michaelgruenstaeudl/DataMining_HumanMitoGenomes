#!/usr/bin/env nextflow

process merge_pacvr_evenness_figures {

    tag "${sample_id}"

    // container 'evolbioinfo/schmutzi:v1.5.6'
    publishDir "results/${sample_id}", mode: 'copy', overwrite: true
    container "texlive/texlive:latest"

    input:
    tuple val(sample_id), val(accession_number), path(pacvr_evenness_initial_figure), path(pacvr_evenness_figure_after_contamination_removal), path(merge_pacvr_evenness_fig_script)

    output:
    // tuple val(sample_id), path("pacvr/${sample_id}_original_CoverageViz.pdf"), emit: evenness_output_ch
    path "combined_coverage_visualization/*"

    script:

    """
    mkdir -p combined_coverage_visualization
    bash ${merge_pacvr_evenness_fig_script} ${sample_id} ${accession_number} ${pacvr_evenness_initial_figure} ${pacvr_evenness_figure_after_contamination_removal}
    """
}
