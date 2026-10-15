# Pipeline usage

The entry point is [`nextflow_workflow/main.nf`](../nextflow_workflow/main.nf).
Run commands from the `nextflow_workflow/` directory so relative default paths
resolve against the project.

## Basic execution

```bash
cd nextflow_workflow
nextflow run main.nf -profile local --outdir results
```

Existing command-line parameter names are retained. Input paths may be
overridden without editing configuration files:

```bash
nextflow run main.nf -profile local \
  --csv_file /path/to/samples.csv \
  --fastq_directory /path/to/fastq \
  --gb_directory /path/to/genbank \
  --fasta_directory /path/to/fasta \
  --reference_fasta /path/to/combined_references.fasta \
  --contamination_db /path/to/kraken_db \
  --config_file /path/to/config.txt \
  --rotation_seed_fasta /path/to/rotate_start_seed.txt
```

## Parameters

| Parameter | Type | Default | Purpose |
| --- | --- | --- | --- |
| `outdir` | path | `results` | Published output directory. |
| `validate_inputs` | boolean | `true` | Validate files, directories, and CSV columns before execution. |
| `csv_file` | path | project-specific | Sample CSV containing `SRA_Run`, `Nucleotide_AccessionID`, and `LibraryLayout`. |
| `fastq_directory` | path | `../genome_data/SRA_data` | Raw FASTQ directory. |
| `gb_directory` | path | `../genome_data/complete_genome/genbank` | GenBank input directory. |
| `fasta_directory` | path | `../genome_data/complete_genome/fasta` | FASTA input directory. |
| `reference_fasta` | path | `reference_genome/combined_references.fasta` | Mapping reference FASTA. |
| `contamination_db` | path | `db` | Kraken2 contamination database. |
| `config_file` | path | `config.txt` | NOVOPlasty configuration. |
| `rotation_seed_fasta` | path | `rotate_start_seed.txt` | Rotation seed sequence. |
| `seed_mito` | path | `Seed_mito.fasta` | NOVOPlasty seed genome. |
| `index_directory` | path | `reference_genome/index` | Reference index location retained for compatibility. |
| `eteDbDir` | path | `myetetoolkit` | ETE toolkit data directory. |
| `calculate_sequence_length_threshold_script` | path | `lib/compute_sequence_length_statistics.py` | Computes read-length cutoffs. |
| `visualize_file_size_changes_script` | path | `lib/visualize_file_size_changes.py` | Generates file-size plots. |
| `genome_comparison_plot_script` | path | `lib/genome_comparison_plot.py` | Generates genome comparison plots. |
| `merge_pacvr_evenness_figure_script` | path | `lib/merge_pacvr_evenness_figures.sh` | Merges coverage figures when evenness is enabled. |
| `trim_mitogenome_duplicate_script` | path | `lib/trim_terminal_suffix_dup.py` | Removes terminal duplicated sequence. |
| `contamination_confidence_threshold` | number | `0.8` | Kraken2 confidence threshold, from 0 through 1. |
| `is_calculate_evenness` | boolean | `true` | Enable initial and final coverage evenness calculations. |
| `rotation_mismatch_threshold` | integer | `5` | Maximum mismatch threshold for genome rotation. |

The authoritative machine-readable parameter contract is
[`nextflow_schema.json`](../nextflow_workflow/nextflow_schema.json).

## Profiles

| Profile | Execution | Use |
| --- | --- | --- |
| `local` | Host executor; containers disabled | Development with the existing environment. |
| `apptainer` | Local executor with Apptainer enabled | Containerized execution where Apptainer is installed. |
| `docker` | Local executor with Docker enabled | Containerized execution where Docker is installed. |
| `test` | Includes `conf/test.config` | Resource-light validation using supplied real inputs. |

Profiles do not install software or create databases. Verify the required
executables, container runtime, reference files, and contamination database
before running.

## Resume

Use `-resume` after correcting a failed task or restarting an interrupted run:

```bash
nextflow run main.nf -profile local -resume --outdir results
```
