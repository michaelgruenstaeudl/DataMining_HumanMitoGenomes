#!/bin/bash

SOURCE_DIR=$1
OUTPUT_CSV=$2

echo "SRA_Number,Topology,Sequence_Length" >"$OUTPUT_CSV"

for dir in "$SOURCE_DIR"/*/; do
    parent=$(basename "$dir")

    if [[ ! "$parent" =~ ^[[:alnum:]]+$ ]]; then
        echo "Skipping $parent"
        continue
    fi

    header=$(sed -n '1p' "$dir"/megahit/mito_denovo.megahit.result/mito_denovo.megahit.mitogenome.fa)
    topology=$(
        echo "$header" |
            awk '
            {
                for(i=1;i<=NF;i++) 
                    if($i ~ /^topology=/)
                        {
                            split($i,a,"="); 
                            print a[2]
                        }
            }
            '
    )

    sequence_length=$(sed -n '2p' "$dir"/megahit/mito_denovo.megahit.result/mito_denovo.megahit.mitogenome.fa | wc -m)
    echo "$parent,$topology,$sequence_length" >>"$OUTPUT_CSV"
done
