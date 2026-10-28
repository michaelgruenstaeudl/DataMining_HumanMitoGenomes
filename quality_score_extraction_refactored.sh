#!/usr/bin/env bash

# Extract initial and post-QC FastQC quality summaries for each SRA run.
# The original script is intentionally left unchanged; this is a drop-in,
# configurable refactor of its behavior.

set -Eeuo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

csv_file="$script_dir/sensitive_data/test_single.csv"
results_dir="$script_dir/nextflow_workflow/results_test"
python_command="${PYTHON:-python3}"
helper_script="$script_dir/nextflow_workflow/lib/calculate_mean_phred_from_fastqc.py"

usage() {
    cat <<'EOF'
Usage: quality_score_extraction_refactored.sh [OPTIONS]

Options:
  --csv FILE       Input CSV (default: sensitive_data/test_single.csv)
  --results DIR    Results directory (default: nextflow_workflow/results_test)
  --python FILE    Python executable (default: PYTHON or python3)
  -h, --help       Show this help message
EOF
}

while (($#)); do
    case "$1" in
        --csv)
            [[ $# -ge 2 ]] || { echo "--csv requires a value" >&2; exit 2; }
            csv_file=$2
            shift 2
            ;;
        --results)
            [[ $# -ge 2 ]] || { echo "--results requires a value" >&2; exit 2; }
            results_dir=$2
            shift 2
            ;;
        --python)
            [[ $# -ge 2 ]] || { echo "--python requires a value" >&2; exit 2; }
            python_command=$2
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

[[ -f "$csv_file" ]] || { echo "Input CSV not found: $csv_file" >&2; exit 1; }
[[ -f "$helper_script" ]] || { echo "Helper script not found: $helper_script" >&2; exit 1; }
command -v "${python_command%% *}" >/dev/null 2>&1 || {
    echo "Python command not found: $python_command" >&2
    exit 1
}

mkdir -p "$results_dir"
output_file="$results_dir/quality_scores.csv"
temporary_output=$(mktemp "${output_file}.tmp.XXXXXX")
cleanup() {
    rm -f -- "$temporary_output"
}
trap cleanup EXIT

printf '%s\n' \
    'RUN,LibraryLayout,Initial_Mean_Phred,Initial_Read_Count,Final_Mean_Phred,Final_Read_Count' \
    >"$temporary_output"

calculate_quality() {
    local run=$1
    local directory=$2
    local -a fastqc_files=()

    if [[ ! -d "$directory" ]]; then
        echo "FastQC directory not found for $run: $directory" >&2
        return 1
    fi

    mapfile -d '' -t fastqc_files < <(
        find "$directory" -maxdepth 1 -type f -name "${run}_*_fastqc.zip" -print0 \
            | sort -z
    )

    if ((${#fastqc_files[@]} == 0)); then
        echo "No FastQC ZIP files found for $run in $directory" >&2
        return 1
    fi
    if ((${#fastqc_files[@]} > 2)); then
        echo "Expected one or two FastQC ZIP files for $run, found ${#fastqc_files[@]}" >&2
        return 1
    fi

    "$python_command" "$helper_script" --fastqc_zip "${fastqc_files[@]}"
}

failures=0
line_number=1
{
    IFS= read -r header || {
        echo "Input CSV is empty: $csv_file" >&2
        exit 1
    }
    [[ -n "$header" ]] || {
        echo "Input CSV has an empty header: $csv_file" >&2
        exit 1
    }

    while IFS=, read -r _accession _sra_uid sra_run library_layout || [[ -n "$sra_run" ]]; do
        ((line_number += 1))
        [[ -n "$sra_run" ]] || continue

        initial_dir="$results_dir/$sra_run/qc/initial_fastqc"
        final_dir="$results_dir/$sra_run/qc/fastqc_after_quality_control"
        initial_result=''
        final_result=''

        echo "Processing $sra_run ($library_layout)"

        if ! initial_result=$(calculate_quality "$sra_run" "$initial_dir"); then
            echo "Initial quality extraction failed for $sra_run (CSV line $line_number)" >&2
            failures=$((failures + 1))
        fi
        if ! final_result=$(calculate_quality "$sra_run" "$final_dir"); then
            echo "Final quality extraction failed for $sra_run (CSV line $line_number)" >&2
            failures=$((failures + 1))
        fi

        IFS=, read -r initial_mean initial_count <<<"$initial_result"
        IFS=, read -r final_mean final_count <<<"$final_result"
        printf '%s,%s,%s,%s,%s,%s\n' \
            "$sra_run" "$library_layout" \
            "${initial_mean:-}" "${initial_count:-}" \
            "${final_mean:-}" "${final_count:-}" \
            >>"$temporary_output"
    done
} <"$csv_file"

mv -- "$temporary_output" "$output_file"
trap - EXIT

if ((failures)); then
    echo "Wrote $output_file with $failures extraction failure(s)." >&2
    exit 1
fi

echo "Wrote $output_file"
