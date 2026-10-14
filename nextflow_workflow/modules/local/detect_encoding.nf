#!/usr/bin/env nextflow

process detectEncoding {

    tag "${meta.id}"

    input:
    tuple val(meta), path(fastqc_report)

    output:
    tuple val(meta), env("encoding_int"), emit: encoding_output_channel

    script:
    """
    ENCODING=\$(unzip -p ${fastqc_report} \
        "*/fastqc_data.txt" | grep "Encoding" | cut -f2)

    if echo "\$ENCODING" | grep -qiE "Illumina 1\\.[3-7]|Solexa"; then
        encoding_int=64
    else
        encoding_int=33
    fi
    echo \$encoding_int
    """
}
