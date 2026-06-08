#!/usr/bin/env nextflow

process extract_seed {

    tag "${sample_id}"

    publishDir "${parent_output_dir}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(input_file)
    val parent_output_dir

    output:
    // val sample_id, emit: sample_id
    tuple val(sample_id), path("seed_output/${sample_id}_seed.fasta", optional: true), emit: seed_output_channel

    script:
    // def merged_file_name = "merged_${sample_id}"
    // # pear -f "${input_file1}" -r "${input_file2}" -o "${merged_file_name}"
    def seed_file = "seed_output/${sample_id}_seed.fasta"
    """
    mkdir seed_output
    
    header=\$(head -n 1 "${input_file}")
    sequence=\$(sed -n '2p' "${input_file}")

    # Clean header name (remove '@' and extra info)
    header_name=\$(echo "\$header" | cut -d' ' -f1 | sed 's/@//')
        echo ">\$header_name" > "${seed_file}"
    echo "\$sequence" >> "${seed_file}"
    """
}
