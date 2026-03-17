#!/usr/bin/env nextflow


process trim_mitogenome_duplicate {
    label 'python_env'
    tag "${sample_id}"
    publishDir "${parent_output_dir}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(genome_fasta), path(script_to_trim_mitogenome_duplicate)
    val parent_output_dir

    output:
    tuple val(sample_id), path("${trimmed_output_fasta}"), emit: trimmed_fasta_ch

    script:
    directory_path = "trimmed_assembled_genomes"
    trimmed_output_fasta = "${directory_path}/${sample_id}_trimmed.fa"
    """
    mkdir -p ${directory_path}
    python ${script_to_trim_mitogenome_duplicate} ${genome_fasta}  > ${trimmed_output_fasta}
    """
}
