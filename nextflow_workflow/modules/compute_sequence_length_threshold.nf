#!/usr/bin/env nextflow

process calculate_sequence_length_threshold {

    tag "${sample_id}"
    container 'community.wave.seqera.io/library/python:3.13.7--b46958bde3c7e023'

    input:
    tuple val(sample_id), path(sra_file_path)

    output:
    tuple val(sample_id), env("lower_cutoff"), env("upper_cutoff"), emit: length_cutoffs

    script:
    """
    cutoffs=\$(python ${workflow.projectDir}/bin/compute_sequence_length_statistics.py --fastq_file ${sra_file_path})
    lower_cutoff=\$(echo \$cutoffs | cut -d ',' -f 1)
    upper_cutoff=\$(echo \$cutoffs | cut -d ',' -f 2)
    echo "\$lower_cutoff \$upper_cutoff"
    """
}
