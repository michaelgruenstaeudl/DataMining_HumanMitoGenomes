#!/usr/bin/env nextflow


process trim_mitogenome_duplicate {
    label 'python_env'
    tag "${meta.id}"

    input:
    tuple val(meta), path(genome_fasta), path(script_to_trim_mitogenome_duplicate)

    output:
    tuple val(meta), path("${trimmed_output_fasta}"), emit: trimmed_fasta_ch

    script:
    directory_path = "trimmed_assembled_genomes"
    trimmed_output_fasta = "${directory_path}/${meta.id}_trimmed.fa"
    """
    set -euo pipefail
    mkdir -p ${directory_path}
    python ${script_to_trim_mitogenome_duplicate} ${genome_fasta}  > ${trimmed_output_fasta}
    """
}
