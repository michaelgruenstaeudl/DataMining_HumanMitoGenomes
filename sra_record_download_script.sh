#!/bin/bash

csv_file="sensitive_data/Sample_of_25_PAIRED_LibraryLayout_SRA_Records_Associated_with_Unique_Mitogenomes_except_y_chromosome.csv"
temp_file=$(mktemp)
# Skip header and write to temp file
tail -n +2 "$csv_file" >"$temp_file"
count=1
while IFS=, read -r RUN _; do
    echo "SRA: $RUN"
    # echo "Sample Name: $Samples"
    # echo "Accession: $Accession"
    fasterq-dump -O "SRA_data" -p -S "$RUN"
    # efetch -db nuccore -id "$Accession" -format fasta >Data/Analysis2/step7_assembly_quality_metric/nucleotide_complete_genome/fasta/"$Accession".fasta
    echo "$count" completed
    count=$((count + 1))
    echo "----------------------"
done <"$temp_file"

# Clean up temporary file after processing
rm "$temp_file"
