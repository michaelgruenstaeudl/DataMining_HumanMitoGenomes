process extract_mapped_fastq {

    label "samtools_env"
    tag "${sample_id}"

    publishDir "${parent_output_dir}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(sam_file)
    val parent_output_dir

    output:
    tuple val(sample_id), path("${sample_id}_mapped_R1.fq"), path("${sample_id}_mapped_R2.fq"), emit: extract_mapped_fastq_output_ch

    script:
    """
    samtools view -b -F 4 ${sam_file} |
    samtools sort -o ${sample_id}.bam

    samtools index ${sample_id}.bam

    samtools fastq \
        -1 ${sample_id}_mapped_R1.fq \
        -2 ${sample_id}_mapped_R2.fq \
        -0 /dev/null \
        -s /dev/null \
        -n \
        ${sample_id}.bam
    """
}
