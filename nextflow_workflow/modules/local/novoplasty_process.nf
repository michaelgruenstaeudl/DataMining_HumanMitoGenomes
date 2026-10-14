#!/usr/bin/env nextflow

process novoplast_process {

    tag "${meta.id}"

    input:
    tuple val(meta), path(input_file_list), val(read_length), val(insert_size), path(seed_mito_path), path(config_file_path)

    output:
    tuple val(meta), path("novoplasty_output/*_mito_extend.fasta"), emit: novoplasty_output_ch
    tuple val(meta), path("novoplasty_output"), emit: novoplasty_output_publish_dir

    script:
    """
    mkdir novoplasty_output

    sed -i "s|^Read Length.*|Read Length            = ${read_length}|" ${config_file_path}
    sed -i "s|^Insert size.*|Insert size            = ${insert_size}|" ${config_file_path}
    sed -i "s|^Seed Input.*|Seed Input            = ${seed_mito_path}|" ${config_file_path}
    sed -i "s|Extend seed directly.*|Extend seed directly  = yes|" ${config_file_path}
    if [ "${meta.single_end}" = "true" ]; then
    sed -i "s|^Single/Paired.*|Single/Paired         = SE|" ${config_file_path}
        sed -i "s|^Combined reads.*|Combined reads        = ${input_file_list[0]}|" ${config_file_path}
    else
        sed -i "s|^Single/Paired.*|Single/Paired         = PE|" ${config_file_path}
        sed -i "s|^Forward reads.*|Forward reads         = ${input_file_list[0]}|" ${config_file_path}
        sed -i "s|^Reverse reads.*|Reverse reads         = ${input_file_list[1]}|" ${config_file_path}
    fi
    
    sed -i "s|^Output path.*|Output path            = novoplasty_output/|" ${config_file_path}

    NOVOPlasty4.3.5.pl -c ${config_file_path}
    """
}
