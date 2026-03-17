#!/usr/bin/env nextflow


process trim_mitogenome_duplicate {
    label 'python_env'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(genome_fasta), val(script_to_trim_mitogenome_duplicate)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fa"), emit: trimmed_fasta_ch
    path "${sample_id}_trim.log"

    script:
    trimmed_output_fasta = "trimmed_assembled_genomes/${sample_id}_trimmed.fa"
    """
    mkdir -p trimmed_assembled_genomes
    python ${script_to_trim_mitogenome_duplicate} ${genome_fasta}  > ${trimmed_output_fasta} 2>&1 ${sample_id}_trim.log
    """
}
