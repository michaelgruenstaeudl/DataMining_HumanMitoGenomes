# Testing and validation

The lightweight test profile is defined in
[`nextflow_workflow/conf/test.config`](../nextflow_workflow/conf/test.config).
It uses one CPU, 2 GB memory, a 30-minute task limit, disables evenness
calculation, and writes to `results_test/`. It still requires a real CSV and
the corresponding input files and software; no synthetic scientific dataset
is assumed.

## Static checks

From the repository root:

```bash
bash tests/test_pipeline.sh
bash -n nextflow_workflow/SCRIPT_collect_genomes.sh
bash -n nextflow_workflow/SCRIPT_collect_genomes_rotation_log.sh
git diff --check
```

`tests/test_pipeline.sh` skips Nextflow checks when Nextflow is not installed.
Set `NEXTFLOW_TEST_CSV` to a representative CSV to enable the stub-run check.

When Nextflow is installed, run the configuration and syntax checks from the
workflow directory:

```bash
cd nextflow_workflow
nextflow -version
nextflow config -profile local
nextflow config -profile test
nextflow run main.nf -profile local -stub-run --outdir results_stub
```

If nf-core tools are installed, run the supported lint command for that
version:

```bash
nf-core lint
```

## Regression run

After the small dataset and environment are available, record the baseline
command and compare the original and refactored runs using the same CSV,
parameters, output file names, sample directories, file counts, and checksums
for important FASTA, BAM, and report files:

```bash
cd nextflow_workflow
nextflow run main.nf -profile test \
  --csv_file /path/to/small.csv \
  --outdir results_test
nextflow run main.nf -profile test -resume \
  --csv_file /path/to/small.csv \
  --outdir results_test
```

The regression is not claimed as passed until those commands complete and the
published outputs have been compared with the baseline.
