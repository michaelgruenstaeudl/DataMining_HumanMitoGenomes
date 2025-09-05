#!/usr/bin/env nextflow

// output:
//     path "novoplasty_outptut/*"
process novoplast_process {

    tag "${sample_id}"

    publishDir "results/${sample_id}", mode: 'copy'
    container 'community.wave.seqera.io/library/novoplasty:4.3.5--d66ab53450fa5022'

    input:
    tuple val(sample_id), path(input_file1), path(input_file2), val(read_length), val(insert_size), path(seed_mito_path), path(config_file_path)

    output:
    path "novoplasty_output/*"

    script:
    """
    mkdir novoplasty_output
    echo read_length: ${read_length}
    echo insert_size: ${insert_size}
    echo input_file1: ${input_file1}
    echo input_file2: ${input_file2}
    echo seed_mito_path: ${seed_mito_path}
    echo config_file_path: ${config_file_path}
    pwd ${config_file_path}
    sed -i "s|^Read Length.*|Read Length            = ${read_length}|" ${config_file_path}
    sed -i "s|^Insert size.*|Insert size            = ${insert_size}|" ${config_file_path}
    sed -i "s|^Seed Input.*|Seed Input            = ${seed_mito_path}|" ${config_file_path}
    sed -i "s|^Forward reads.*|Forward reads            = ${input_file1}|" ${config_file_path}
    sed -i "s|^Reverse reads.*|Reverse reads            = ${input_file2}|" ${config_file_path}
    sed -i "s|^Output path.*|Output path            = novoplasty_output/|" ${config_file_path}

    NOVOPlasty4.3.5.pl -c ${config_file_path}
    """
}
