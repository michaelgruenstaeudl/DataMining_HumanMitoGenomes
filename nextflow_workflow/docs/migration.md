# Migration notes

This repository uses an incremental nf-core-inspired refactor. It is not a
template-generated nf-core pipeline and does not add external nf-core module
dependencies.

## Changes by stage

1. **Configuration and validation**: configuration was split into `conf/`, input
   validation was added to `main.nf`, and the existing parameter names were
   retained.
2. **Module layout**: active processes were moved to
   `nextflow_workflow/modules/local/`, with logical workflow groups in
   `nextflow_workflow/subworkflows/local/`.
3. **Input contract**: sample-specific channels use `tuple(meta, ...)`, with
   stable `meta.id`, accession, and single-end metadata where available.
4. **Subworkflow boundaries**: input preparation, QC, contamination and
   mapping, assembly and normalization, rotation, haplotype, and reporting
   responsibilities are exposed through explicit subworkflows.
5. **Error handling**: command failures now fail their tasks, and required
   output checks report missing QC reads, alignments, or assemblies.
6. **Output publication**: output destinations are grouped under per-sample
   stage directories, with aggregate reports and execution logs at the top
   level.

## Compatibility

The CSV format, parameter names, sample identifiers, and biological command
implementations are preserved. Existing scripts that consumed rotated genomes
were updated for the new rotation directory. No compatibility alias for the
old flat output layout is created because repository searches found no active
notebook or script that requires it.

The current small-dataset regression comparison remains pending until the
dataset and a working Nextflow environment are available.
