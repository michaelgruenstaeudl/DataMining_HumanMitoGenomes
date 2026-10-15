process extract_assembled_genome_metadata_info {
    input:
    tuple val(meta), path(assembled_genome_fasta)

    output:
    tuple val(meta), env("topology"), env("sequence_length"), emit: metadata_info_out_ch

    script:
    """
    set -euo pipefail
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

    sequence_length=\$(awk '!/^>/ {gsub(/[[:space:]]/, ""); total += length} END {print total + 0}' ${assembled_genome_fasta})

    if [[ "\$sequence_length" -le 0 ]]; then
        echo "No valid mitochondrial sequence was produced for ${meta.id}." >&2
        exit 1
    fi

    echo "\$topology \$sequence_length"
    """
}
