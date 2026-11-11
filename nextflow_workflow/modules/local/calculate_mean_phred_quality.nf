#!/usr/bin/env nextflow

process calculate_mean_phred_quality {

    label 'python_env'
    tag "${meta.id}"

    input:
    tuple val(meta), path(fastqc_zip), path(calculate_mean_phred_quality_script)

    output:
    tuple val(meta), env("phred_quality_score"), env("total_reads"), emit: phred_quality

    script:
    """
    set -euo pipefail
    quality_result=\$(python ${calculate_mean_phred_quality_script} --fastqc_zip ${fastqc_zip})
    IFS=',' read -r phred_quality_score total_reads <<< "\$quality_result"

    if [[ ! "\$phred_quality_score" =~ ^[0-9]+([.][0-9]+)?\$ ]]; then
        echo "Invalid mean Phred quality returned for ${meta.id}: \$quality_result" >&2
        exit 1
    fi
    if [[ ! "\$total_reads" =~ ^[1-9][0-9]*\$ ]]; then
        echo "Invalid read count returned for ${meta.id}: \$quality_result" >&2
        exit 1
    fi

    export phred_quality_score total_reads
    printf '%s,%s\\n' "\$phred_quality_score" "\$total_reads"
    """
}
