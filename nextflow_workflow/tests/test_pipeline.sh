#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
workflow_dir="${project_root}/nextflow_workflow"

if ! command -v nextflow >/dev/null 2>&1; then
    echo "SKIP: nextflow is not installed"
    exit 0
fi

cd "${workflow_dir}"
nextflow config -profile local >/dev/null
nextflow config -profile test >/dev/null
echo "Nextflow configuration checks: OK"

if [[ -n "${NEXTFLOW_TEST_CSV:-}" ]]; then
    nextflow run main.nf -profile test -stub-run \
        --csv_file "${NEXTFLOW_TEST_CSV}" \
        --outdir "${NEXTFLOW_TEST_OUTDIR:-results_stub}"
else
    echo "SKIP: set NEXTFLOW_TEST_CSV to run the stub regression check"
fi
