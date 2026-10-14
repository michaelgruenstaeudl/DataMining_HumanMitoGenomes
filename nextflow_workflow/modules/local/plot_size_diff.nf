#!/usr/bin/env nextflow

process plot_size_diff {
    label 'python_env'

    input:
    tuple path(csv_file_path), path(directory_name), path(visualize_file_size_changes_script_ch)

    script:
    """
    python ${visualize_file_size_changes_script_ch} --csv_file ${csv_file_path} --directory_name ${directory_name}
    """
}
