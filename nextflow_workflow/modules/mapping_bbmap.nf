#!/usr/bin/env nextflow

process mapping_bbmap {

    label "bbmap_env"
    tag "${sample_id}"

    // container 'fastq_sifter'

    publishDir "${parent_output_dir}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(input_file1), path(input_file2)
    path index_dir
    val parent_output_dir

    output:
    tuple val(sample_id), path("bbmap_mapped_output/${sample_id}_mapped_R1.fq"), path("bbmap_mapped_output/${sample_id}_mapped_R2.fq"), emit: mapping_process_output
    path "bbmap_mapped_output/${sample_id}_unmapped_R1.fq"
    path "bbmap_mapped_output/${sample_id}_unmapped_R2.fq"

    script:
    def avail_mem = (task.memory.toGiga() * 0.85).intValue()
    """
    mkdir bbmap_mapped_output
    bbmap.sh \
        -Xmx${avail_mem}g \
        path=${index_dir} \
        in1=${input_file1} in2=${input_file2} \
        outm1=bbmap_mapped_output/${sample_id}_mapped_R1.fq outm2=bbmap_mapped_output/${sample_id}_mapped_R2.fq \
        outu1=bbmap_mapped_output/${sample_id}_unmapped_R1.fq outu2=bbmap_mapped_output/${sample_id}_unmapped_R2.fq
    """
}
