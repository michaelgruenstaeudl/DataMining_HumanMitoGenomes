process extract_haplotype_info {
    input:
    tuple val(meta), path(haplotype_csv)

    output:
    // tuple val(meta), env("topology"), env("sequence_length"), emit: metadata_info_out_ch
    tuple val(meta), env("haplotype"), env("quality"), emit: haplotype_info_out_ch

    script:
    """
    echo "Starting haplotype information extraction for sample: ${meta.id}"
    # cat ${haplotype_csv}
    # Extracting haplotype information from the CSV file
    haplotype=\$(awk -F' ' 'NR==2 {print \$2}' ${haplotype_csv})
    quality=\$(awk -F' ' 'NR==2 {print \$4}' ${haplotype_csv})
    echo \$haplotype \$quality
    """
}
