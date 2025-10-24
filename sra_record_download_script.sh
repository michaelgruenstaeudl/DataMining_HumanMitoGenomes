#!/bin/bash

csv_file="sensitive_data/sample_nucleotide_sra_mapped_infor_with_run.csv"
temp_file=$(mktemp)
# Skip header and write to temp file
tail -n +2 "$csv_file" >"$temp_file"
count=1
while IFS=, read -r Accession _SRA_Uid RUN _LibraryLayout; do
    echo "SRA: $RUN"
    # echo "Sample Name: $Samples"
    echo "Accession: $Accession"
    fasterq-dump -O "genome_data/SRA_data" -p -S "$RUN"
    efetch -db nuccore -id "$Accession" -format fasta >genome_data/complete_genome/fasta/"$Accession".fasta
    efetch -db nuccore -id "$Accession" -format gb >genome_data/complete_genome/genbank/"$Accession".gb

    echo "$count" completed
    count=$((count + 1))
    echo "----------------------"
done <"$temp_file"

# Clean up temporary file after processing
rm "$temp_file"
