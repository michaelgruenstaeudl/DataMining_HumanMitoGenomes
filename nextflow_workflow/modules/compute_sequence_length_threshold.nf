#!/usr/bin/env nextflow

process calculate_sequence_length_threshold {

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(sra_file_path)

    output:
    tuple val(sample_id), env("lower_cutoff"), env("upper_cutoff"), emit: length_cutoffs

    script:
    """
    cutoffs=\$(python3 /home/b_thapamagar/BioInformatics/NCBIrecordMining/nextflow_workflow/modules/compute_sequence_length_statistics.py --fastq_file ${sra_file_path})
    lower_cutoff=\$(echo \$cutoffs | cut -d ',' -f 1)
    upper_cutoff=\$(echo \$cutoffs | cut -d ',' -f 2)
    echo "\$lower_cutoff \$upper_cutoff"
    """
}
