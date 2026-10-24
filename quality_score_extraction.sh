#!/bin/bash
csv_file="sensitive_data/test_single.csv"
temp_file=$(mktemp)
results_dir="nextflow_workflow/results_test"

tail -n +2 "$csv_file" >"$temp_file"

echo "RUN,LibraryLayout,Initial_Mean_Phred,Initial_Read_Count,Final_Mean_Phred,Final_Read_Count" >"${results_dir}/quality_scores.csv"
while IFS=, read -r _Accession _SRA_Uid SRA_RUN LibraryLayout; do
    initial_fastqc_dir="${results_dir}/${SRA_RUN}/qc/initial_fastqc"
    final_fastqc_dir="${results_dir}/${SRA_RUN}/qc/fastqc_after_quality_control"
    echo "RUN: $SRA_RUN"
    echo "Initial FastQC dir: $initial_fastqc_dir"
    echo "Final FastQC dir: $final_fastqc_dir"
    mapfile -t initial_fastqc_files < <(
        find "$initial_fastqc_dir" -maxdepth 1 -type f \
            -name "${SRA_RUN}"_*_fastqc.zip
    )

    echo "${initial_fastqc_files[@]}"

    # if ((${#fastqc_files[@]} == 2)); then
    if ((${#initial_fastqc_files[@]} >= 1)); then
        initial_result=$(python nextflow_workflow/lib/calculate_mean_phred_from_fastqc.py \
            --fastqc_zip "${initial_fastqc_files[@]}")
        echo "Initial result: $initial_result"
    else
        echo "Expected 1 or 2 FastQC ZIP files for ${SRA_RUN}, found ${#initial_fastqc_files[@]}" >&2
    fi

    mapfile -t final_fastqc_files < <(
        find "$final_fastqc_dir" -maxdepth 1 -type f \
            -name "${SRA_RUN}"_*_fastqc.zip
    )
    echo "${final_fastqc_files[@]}"
    if ((${#final_fastqc_files[@]} >= 1)); then
        final_result=$(python nextflow_workflow/lib/calculate_mean_phred_from_fastqc.py \
            --fastqc_zip "${final_fastqc_files[@]}")
        echo "Final result: $final_result"
    else
        echo "Expected 1 or 2 FastQC ZIP files for ${SRA_RUN}, found ${#final_fastqc_files[@]}" >&2
    fi

    echo "${SRA_RUN},${LibraryLayout},${initial_result},${final_result}" >>"${results_dir}/quality_scores.csv"

    echo "-----------------------------"
done <"$temp_file"

# Clean up temporary file after processing
rm "$temp_file"
