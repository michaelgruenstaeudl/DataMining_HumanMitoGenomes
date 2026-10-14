#!/usr/bin/env nextflow

process mapping_bbmap {

    label "bbmap_env"
    tag "${meta.id}"

    input:
    tuple val(meta), path(input_file_list)
    path index_dir

    output:
    tuple val(meta), path("bbmap_mapped_output/${meta.id}_mapped*.fq"), emit: mapping_process_output
    tuple val(meta), path("bbmap_mapped_output"), optional: true, emit: mapping_process_publish_dir

    script:
    def avail_mem = (task.memory.toGiga() * 0.85).intValue()
    def input_file = meta.single_end
        ? "in=${input_file_list[0]}"
        : "in=${input_file_list[0]} in2=${input_file_list[1]}"
    def mapped_output = meta.single_end
        ? "outm=bbmap_mapped_output/${meta.id}_mapped.fq"
        : "outm1=bbmap_mapped_output/${meta.id}_mapped_R1.fq outm2=bbmap_mapped_output/${meta.id}_mapped_R2.fq"
    def unmapped_output = meta.single_end
        ? "outu=bbmap_mapped_output/${meta.id}_unmapped.fq"
        : "outu1=bbmap_mapped_output/${meta.id}_unmapped_R1.fq outu2=bbmap_mapped_output/${meta.id}_unmapped_R2.fq"
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
