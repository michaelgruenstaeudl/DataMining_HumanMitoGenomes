#!/usr/bin/env nextflow

process novoplast_process {

    tag "${sample_id}"

    // container 'community.wave.seqera.io/library/novoplasty:4.3.5--d66ab53450fa5022'
    publishDir "${parent_output_dir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(input_file_list), val(is_single_end), val(read_length), val(insert_size), path(seed_mito_path), path(config_file_path)
    val parent_output_dir

    output:
    path "novoplasty_output/*"

    script:
    """
    mkdir novoplasty_output

    sed -i "s|^Read Length.*|Read Length            = ${read_length}|" ${config_file_path}
    sed -i "s|^Insert size.*|Insert size            = ${insert_size}|" ${config_file_path}
    sed -i "s|^Seed Input.*|Seed Input            = ${seed_mito_path}|" ${config_file_path}
    sed -i "s|Extend seed directly.*|Extend seed directly  = yes|" ${config_file_path}
    if [ "${is_single_end}" = "true" ]; then
        sed -i "s|^Single/Paired.*|Single/Paired         = SE|" ${config_file_path}
        sed -i "s|^Combined reads.*|Combined reads        = ${input_file_list[1]}|" ${config_file_path}
    else
        sed -i "s|^Single/Paired.*|Single/Paired         = SE|" ${config_file_path}
        sed -i "s|^Forward reads.*|Forward reads            = ${input_file_list[0]}|" ${config_file_path}
        sed -i "s|^Reverse reads.*|Reverse reads            = ${input_file_list[1]}|" ${config_file_path}
    fi
    
    sed -i "s|^Output path.*|Output path            = novoplasty_output/|" ${config_file_path}

    NOVOPlasty4.3.5.pl -c ${config_file_path}
    """
}
