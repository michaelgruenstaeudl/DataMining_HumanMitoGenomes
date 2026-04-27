process extract_haplotype_info {

    publishDir "${parent_output_dir}/reports", mode: 'copy'

    input:
    tuple val(sample_id), path(haplotype_csv)
    val parent_output_dir

    output:
    // tuple val(sample_id), env("topology"), env("sequence_length"), emit: metadata_info_out_ch
    tuple val(sample_id), env("haplotype"), env("quality"), emit: haplotype_info_out_ch

    script:
    """
    echo "Starting haplotype information extraction for sample: ${sample_id}"
    # cat ${haplotype_csv}
    # Extracting haplotype information from the CSV file
    haplotype=\$(awk -F' ' 'NR==2 {print \$2}' ${haplotype_csv})
    quality=\$(awk -F' ' 'NR==2 {print \$4}' ${haplotype_csv})
    echo \$haplotype \$quality
    """
}
