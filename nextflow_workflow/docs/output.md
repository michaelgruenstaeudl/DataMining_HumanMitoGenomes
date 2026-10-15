# Output directory

The default output directory is `results/`. The output declarations are in
[`nextflow_workflow/main.nf`](../nextflow_workflow/main.nf); process work files
remain in Nextflow's work directory.

```text
results/
├── <sample>/
│   ├── qc/
│   ├── contamination/
│   ├── mapping/
│   ├── assembly/
│   ├── normalization/
│   ├── rotation/
│   ├── haplotype/
│   └── reports/
├── reports/
└── logs/
```

## Published content

- `qc/`: FastQC, repaired reads, and quality-control output directories.
- `contamination/`: Kraken2 contamination-removal results.
- `mapping/`: BBMap, SAM, BAM, reference-index, and evenness outputs.
- `assembly/`: MitoZ assembly output.
- `normalization/`: NOVOPlasty and normalized FASTA output.
- `rotation/`: assembled and official rotated genomes, including rotation logs.
- `haplotype/`: sample-level HaploGreP output.
- `<sample>/reports/`: sample-level merged coverage figures.
- `reports/`: cutoff, file-size, assembly, normalization, rotation, and haplotype CSV reports.
- `logs/`: Nextflow execution trace.

File names and sample identifiers remain unchanged. The rotated FASTA path
changed from `<sample>/rotated_output/` to
`<sample>/rotation/rotated_output/`. The repository collection scripts were
updated for this path.
