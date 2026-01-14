process extract_assembled_genome_metadata_info {

    publishDir "${parent_output_dir}/reports", mode: 'copy'

    input:
    tuple val(sample_id), path(assembled_genome_fasta)
    val parent_output_dir

    output:
    tuple val(sample_id), env("topology"), env("sequence_length"), emit: metadata_info_out_ch

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

    sequence_length=\$(sed -n '2p' ${assembled_genome_fasta} | wc -m)

    echo "\$topology \$sequence_length"
    """
}
