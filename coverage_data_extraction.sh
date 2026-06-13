#!/bin/bash
output_dir="coverage_data"
mkdir -p "$output_dir"
initial_phase_output_dir="$output_dir/initial_phase_output"
final_phase_output_dir="$output_dir/final_phase_output"
logs_output_dir="$output_dir/logs"
mkdir -p "$initial_phase_output_dir"
mkdir -p "$final_phase_output_dir"
mkdir -p "$logs_output_dir"

script_name=$(basename "$0" .sh)
run_log_file="$logs_output_dir/${script_name}_$(date +%Y%m%d_%H%M%S).log"
touch "$run_log_file"
exec > >(tee -a "$run_log_file") 2>&1

echo "Logging to $run_log_file"
# bedtools genomecov -ibam sample.bam -d >coverage.txt

for dir in nextflow_workflow/results/*; do
    if [ -d "$dir" ]; then
        sra_name=$(basename "$dir")
        echo "Processing directory: $sra_name"
        if [ -f "$dir/bam_output_initial_phase/$sra_name.bam" ]; then
            initial_coverage_file_name="$sra_name"_initial_coverage.csv
            echo "initial BAM file exists for $sra_name"
            bedtools genomecov -ibam "$dir/bam_output_initial_phase/$sra_name.bam" -d >"$initial_phase_output_dir/$initial_coverage_file_name"
        else
            echo "initial BAM file does not exist for $sra_name"
        fi
        if [ -f "$dir/bam_output_mapping_phase/$sra_name.bam" ]; then
            final_coverage_file_name="$sra_name"_final_coverage.csv
            echo "mapping BAM file exists for $sra_name"
            bedtools genomecov -ibam "$dir/bam_output_mapping_phase/$sra_name.bam" -d >"$final_phase_output_dir/$final_coverage_file_name"

        else
            echo "mapping BAM file does not exist for $sra_name"
        fi
    fi
done
