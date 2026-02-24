#!/bin/bash
# csv_file="sensitive_data/test_single.csv"
csv_file="sensitive_data/SINGLE_layout_sample_nucleotide_sra_mapped_info_with_run.csv"
temp_file=$(mktemp)
output_dir="single_layout_genome_data"
SRA_data_dir="$output_dir/SRA_data"
complete_genome_fasta_dir="$output_dir/complete_genome/fasta"
complete_genome_genbank_dir="$output_dir/complete_genome/genbank"

mkdir -p "$SRA_data_dir"
mkdir -p "$complete_genome_fasta_dir"
mkdir -p "$complete_genome_genbank_dir"
# Skip header and write to temp file
tail -n +2 "$csv_file" >"$temp_file"
count=1
while IFS=, read -r Accession _SRA_Uid RUN LibraryLayout; do
    echo "SRA: $RUN"
    # echo "Sample Name: $Samples"
    echo "Accession: $Accession"

    FASTQ1="${SRA_data_dir}/${RUN}_1.fastq"
    FASTQ2="${SRA_data_dir}/${RUN}_2.fastq"
    FASTA="${complete_genome_fasta_dir}/${Accession}.fasta"
    GB="${complete_genome_genbank_dir}/${Accession}.gb"

    # Download FASTQ if missing
    if [[ "$LibraryLayout" == "SINGLE" ]]; then
        existing=("${SRA_data_dir}/${RUN}"*.fastq)
        if [[ ! -e "${existing[0]}" ]]; then
            echo "Downloading single-end FASTQ for $RUN ..."
            fasterq-dump -O "${SRA_data_dir}" -p -S "$RUN"
        else
            echo "Single-end FASTQ for $RUN already exists, skipping."
        fi
    else
        # Paired-end: check both FASTQ1 and FASTQ2
        if [[ ! -f "$FASTQ1" || ! -f "$FASTQ2" ]]; then
            echo "Downloading paired-end FASTQ for $RUN ..."
            fasterq-dump -O "${SRA_data_dir}" -p -S "$RUN"
        else
            echo "Paired-end FASTQ for $RUN already exists, skipping."
        fi
    fi

    # fasterq-dump -O "${SRA_data_dir}" -p -S "$RUN"

    # Download GenBank if missing
    if [[ ! -f "$FASTA" ]]; then
        echo "Downloading FASTA for $Accession ..."
        efetch -db nuccore -id "$Accession" -format fasta >"$FASTA"
    else
        echo "FASTA for $Accession already exists, skipping."
    fi

    # efetch -db nuccore -id "$Accession" -format fasta >${complete_genome_fasta_dir}/"$Accession".fasta

    # Download GenBank if missing
    if [[ ! -f "$GB" ]]; then
        echo "Downloading GenBank for $Accession ..."
        efetch -db nuccore -id "$Accession" -format gb >"$GB"
    else
        echo "GenBank for $Accession already exists, skipping."
    fi

    # efetch -db nuccore -id "$Accession" -format gb >${complete_genome_genbank_dir}/"$Accession".gb

    echo "$count" completed
    count=$((count + 1))
    echo "----------------------"
done <"$temp_file"

# Clean up temporary file after processing
rm "$temp_file"
