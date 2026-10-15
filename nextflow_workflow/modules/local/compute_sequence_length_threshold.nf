#!/usr/bin/env nextflow

process calculate_sequence_length_threshold {

    label 'python_env'
    // container 'community.wave.seqera.io/library/python:3.13.7--b46958bde3c7e023'
    tag "${meta.id}"

    input:
    tuple val(meta), path(sra_file_path), path(calculate_sequence_length_threshold_script_ch)

    output:
    tuple val(meta), env("lower_cutoff"), env("upper_cutoff"), emit: length_cutoffs

    script:
    """
    set -euo pipefail
    cutoffs=\$(python ${calculate_sequence_length_threshold_script_ch} --fastq_file ${sra_file_path})
    lower_cutoff=\$(echo \$cutoffs | cut -d ',' -f 1)
    upper_cutoff=\$(echo \$cutoffs | cut -d ',' -f 2)
    if [[ -z "\$lower_cutoff" || -z "\$upper_cutoff" ]]; then
        echo "No valid sequence-length cutoffs were produced for ${meta.id}." >&2
        exit 1
    fi
    export lower_cutoff upper_cutoff
    echo "\$lower_cutoff \$upper_cutoff"
    """
}
