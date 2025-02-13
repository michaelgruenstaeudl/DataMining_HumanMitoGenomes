#!/bin/sh

SEARCH_QUERY="\"Homo sapiens\"[ORGN] AND \"complete genome\"[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] \
            NOT (unverified OR \"Homo sp. Altai\" OR \"Denisova hominin\" OR neanderthalensis OR heidelbergensis OR consensus)"

echo "SRA records related using BioSample is being fetched:"

# Define the CSV file path
CSV_FILE1="DATA_SRA_Metadata_linked_through_BioSample.csv"

# Check if the file exists, if not, create it and add headers
if [ -f "$CSV_FILE1" ]; then
    rm "$CSV_FILE1"
fi

echo "PRIMARY_ID,BioSample,BioProject,SRA_ID" >"$CSV_FILE1"

esearch -db nucleotide -query "$SEARCH_QUERY" |
    elink -target biosample |
    elink -target sra |
    efetch -format native -mode xml | #>"$FILE_NAME1"
    xtract \
        -pattern EXPERIMENT_PACKAGE \
        -block EXPERIMENT/IDENTIFIERS/PRIMARY_ID -element PRIMARY_ID \
        -block SAMPLE/IDENTIFIERS/EXTERNAL_ID -if "@namespace" -equals "BioSample" -element EXTERNAL_ID \
        -block STUDY/IDENTIFIERS/EXTERNAL_ID -if "@namespace" -equals "BioProject" -element EXTERNAL_ID \
        -block SUBMISSION/IDENTIFIERS/PRIMARY_ID -element PRIMARY_ID |
    tr '\t' ',' >>"$CSV_FILE1"

# curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&query_key=$QUERY_KEY&WebEnv=$WEBENV&rettype=native&retmode=xml" \
#     >>"$FILE_NAME"
echo "SRA metadata fetched successfully and stored in $CSV_FILE1 file"

# SRA data Fetching from Nucleotide -> BioProject -> SRA database link

echo "SRA records related to BioProject is being fetched:"

# Define the CSV file path
CSV_FILE2="DATA_SRA_Metadata_linked_through_BioProject.csv"

# Check if the file exists, if not, create it and add headers
if [ -f "$CSV_FILE2" ]; then
    rm "$CSV_FILE2"
fi

echo "PRIMARY_ID,BioSample,BioProject,SRA_ID" >"$CSV_FILE2"

esearch -db nucleotide -query "$SEARCH_QUERY" |
    elink -target bioproject |
    elink -target sra |
    efetch -format native -mode xml | #>"$FILE_NAME1"
    xtract \
        -pattern EXPERIMENT_PACKAGE \
        -block EXPERIMENT/IDENTIFIERS/PRIMARY_ID -element PRIMARY_ID \
        -block SAMPLE/IDENTIFIERS/EXTERNAL_ID -if "@namespace" -equals "BioSample" -element EXTERNAL_ID \
        -block STUDY/IDENTIFIERS/EXTERNAL_ID -if "@namespace" -equals "BioProject" -element EXTERNAL_ID \
        -block SUBMISSION/IDENTIFIERS/PRIMARY_ID -element PRIMARY_ID |
    tr '\t' ',' >>"$CSV_FILE2"

echo "SRA metadata fetched successfully and stored in $CSV_FILE2 file"
