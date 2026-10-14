process extract_assembled_genome_metadata_info {
    input:
    tuple val(meta), path(assembled_genome_fasta)

    output:
    tuple val(meta), env("topology"), env("sequence_length"), emit: metadata_info_out_ch

    script:
    """
    header=\$(sed -n '1p' ${assembled_genome_fasta})
    
    topology=\$(
        echo "\$header" |
            awk '
            {
                for(i=1;i<=NF;i++) 
                    if(\$i ~ /^topology=/)
                        {
                            split(\$i,a,"="); 
                            print a[2]
                        }
            }
            '
    )

    sequence_length=\$(grep -v '^>' ${assembled_genome_fasta} | tr -d '[:space:]' | wc -c)

    echo "\$topology \$sequence_length"
    """
}
