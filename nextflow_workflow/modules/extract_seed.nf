#!/usr/bin/env nextflow

process seedExtractionProcess {

    tag "${sample_id}"

    publishDir "results/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(input_file1), path(input_file2)

    output:
    // val sample_id, emit: sample_id
    tuple val(sample_id), path("seed_output/${sample_id}_seed.fasta", optional: true), emit: seed_output_channel

    script:
    def merged_file_name = "merged_${sample_id}"
    def seed_file = "seed_output/${sample_id}_seed.fasta"
    """
    mkdir seed_output
    pear -f "${input_file1}" -r "${input_file2}" -o "${merged_file_name}"
    
    header=\$(head -n 1 "${merged_file_name}.assembled.fastq")
    sequence=\$(sed -n '2p' "${merged_file_name}.assembled.fastq")

    # Clean header name (remove '@' and extra info)
    header_name=\$(echo "\$header" | cut -d' ' -f1 | sed 's/@//')
        echo ">\$header_name" > "${seed_file}"
    echo "\$sequence" >> "${seed_file}"
    """
}

workflow {
    input_ch = Channel.fromFilePairs("/home/b_thapamagar/BioInformatics/DataMining_HumanMitoGenomes/temp_data2/temp" + '/*_{1,2}.fastq')
        .map { sample_id, read -> tuple(sample_id, read[0]) }
        .view { "Input to repair reads: ${it}" }

    seedExtractionProcess(input_ch)
    seedExtractionProcess.out.seed_output_channel.view { "Extracted seed file: ${it}" }
}
