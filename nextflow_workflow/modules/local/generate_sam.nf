#!/usr/bin/env nextflow

process generate_sam {

    tag "${meta.id}"

    input:
    tuple val(meta), path(input_fastq_list), path(fasta_file)
    val phase

    output:
    tuple val(meta), env("exit_code"), emit: sam_status_report
    tuple val(meta), path("sam_output_${phase}/${meta.id}.sam"), emit: sam_output_ch
    tuple val(meta), path("sam_output_${phase}"), emit: sam_output_publish_dir

    script:
    def output_sam_filename = "sam_output_${phase}/${meta.id}.sam"
    input_param = meta.single_end
        ? "-U ${input_fastq_list[0]}"
        : "-1 ${input_fastq_list[0]} -2 ${input_fastq_list[1]}"
    // def output_bam_filename = "bam_output/${sample_id}.bam"
    // #samtools view -@ 8 -Sb -T ${fasta_file} ${output_sam_filename} >${output_bam_filename}

    """
    ## MAPPING ORIGINAL READS
    set +e
    mkdir -p sam_output_${phase}
    
    mkdir -p temp_data2/temp/db
    bowtie2-build ${fasta_file} temp_data2/temp/db/${meta.id}
    bowtie2 -x temp_data2/temp/db/${meta.id} ${input_param} -S ${output_sam_filename} -p 8

    exit_code=\$(echo \$? )
    echo \$exit_code > .exitcode
    echo \$exit_code
    """
}
