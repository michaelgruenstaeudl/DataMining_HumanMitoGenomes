#!/usr/bin/env nextflow

process detectEncoding {

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fastqc_report)

    output:
    // val sample_id, emit: sample_id
    tuple val(sample_id), env("encoding_int"), emit: encoding_output_channel

    script:
    """
    encoding_int=\$(grep -i'encoding' ${fastqc_report} \
    | awk -F':' '{print \$2}' \
    | grep -q "33" && echo "33" || echo "64" )
    echo \$encoding_int
    """
}

workflow {
    input_ch = Channel.of(tuple('sample1', '/home/b_thapamagar/BioInformatics/DataMining_HumanMitoGenomes/temp_data2/temp/ERR322856_1.fastq_trimming_report.txt'))
        .view { "Input to detect encoding: ${it}" }
    detectEncoding(input_ch)
}
