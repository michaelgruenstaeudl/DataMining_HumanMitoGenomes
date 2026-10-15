#!/usr/bin/env nextflow

process merge_pacvr_evenness_figures {

    tag "${meta.id}"

    input:
    tuple val(meta), path(pacvr_evenness_initial_figure), path(pacvr_evenness_figure_after_contamination_removal), path(merge_pacvr_evenness_fig_script)

    output:
    // tuple val(meta), path("pacvr/${meta.id}_original_CoverageViz.pdf"), emit: evenness_output_ch
    tuple val(meta), path("combined_coverage_visualization"), emit: merged_pacvr_figures_output_ch
    tuple val(meta), path("combined_coverage_visualization/*"), emit: merged_pacvr_figures_files

    script:

    """
    set -euo pipefail
    mkdir -p combined_coverage_visualization
    bash ${merge_pacvr_evenness_fig_script} ${meta.id} ${meta.accession} ${pacvr_evenness_initial_figure} ${pacvr_evenness_figure_after_contamination_removal}
    """
}
