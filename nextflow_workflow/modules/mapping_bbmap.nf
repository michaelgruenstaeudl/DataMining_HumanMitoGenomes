#!/usr/bin/env nextflow

process mapping_bbmap {

    label "bbmap_env"
    tag "${sample_id}"

    // container 'fastq_sifter'

    publishDir "${parent_output_dir}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(input_file_list), val(is_single_end)
    path index_dir
    val parent_output_dir

    output:
    tuple val(sample_id), path("bbmap_mapped_output/${sample_id}_mapped*.fq"), val(is_single_end), emit: mapping_process_output
    path "bbmap_mapped_output/${sample_id}_unmapped*.fq", optional: true

    script:
    def avail_mem = (task.memory.toGiga() * 0.85).intValue()
    def input_file = is_single_end
        ? "in=${input_file_list[0]}"
        : "in=${input_file_list[0]} in2=${input_file_list[1]}"
    def mapped_output = is_single_end
        ? "outm=bbmap_mapped_output/${sample_id}_mapped.fq"
        : "outm1=bbmap_mapped_output/${sample_id}_mapped_R1.fq outm2=bbmap_mapped_output/${sample_id}_mapped_R2.fq"
    def unmapped_output = is_single_end
        ? "outu=bbmap_mapped_output/${sample_id}_unmapped.fq"
        : "outu1=bbmap_mapped_output/${sample_id}_unmapped_R1.fq outu2=bbmap_mapped_output/${sample_id}_unmapped_R2.fq"
    """
    mkdir bbmap_mapped_output
    bbmap.sh \
        -Xmx${avail_mem}g \
        path=${index_dir} \
        ${input_file} \
        ${mapped_output} \
        ${unmapped_output}
    """
}
