#!/usr/bin/env nextflow

process detectEncoding {

    tag "${meta.id}"

    input:
    tuple val(meta), path(fastqc_report)

    output:
    tuple val(meta), env("encoding_int"), emit: encoding_output_channel

    script:
    """
    set -euo pipefail
    ENCODING=\$(unzip -p ${fastqc_report} \
        "*/fastqc_data.txt" | awk -F '\\t' '\$1 == "Encoding" { print \$2; exit }')

    if [[ -z "\$ENCODING" ]]; then
        echo "FastQC encoding could not be determined for ${meta.id}." >&2
        exit 1
    fi

    if echo "\$ENCODING" | grep -qiE "Illumina 1\\.[3-7]|Solexa"; then
        encoding_int=64
    else
        encoding_int=33
    fi
    export encoding_int
    echo \$encoding_int
    """
}
