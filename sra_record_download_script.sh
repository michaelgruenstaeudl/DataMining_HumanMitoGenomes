#!/bin/bash

csv_file="sensitive_data/PRJEB4417_nucleotide_sra_mapped_infor_with_run.csv"
temp_file=$(mktemp)
output_dir="PRJEB4417_genome_data"
SRA_data_dir="$output_dir/SRA_data"
complete_genome_fasta_dir="$output_dir/complete_genome/fasta"
complete_genome_genbank_dir="$output_dir/complete_genome/genbank"

mkdir -p "$SRA_data_dir"
mkdir -p "$complete_genome_fasta_dir"
mkdir -p "$complete_genome_genbank_dir"
# Skip header and write to temp file
tail -n +2 "$csv_file" >"$temp_file"
count=1
while IFS=, read -r Accession _SRA_Uid RUN _LibraryLayout; do
    echo "SRA: $RUN"
    # echo "Sample Name: $Samples"
    echo "Accession: $Accession"
    fasterq-dump -O "${SRA_data_dir}" -p -S "$RUN"
    efetch -db nuccore -id "$Accession" -format fasta >${complete_genome_fasta_dir}/"$Accession".fasta
    efetch -db nuccore -id "$Accession" -format gb >${complete_genome_genbank_dir}/"$Accession".gb

    echo "$count" completed
    count=$((count + 1))
    echo "----------------------"
done <"$temp_file"

# Clean up temporary file after processing
rm "$temp_file"
