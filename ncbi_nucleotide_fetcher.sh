#!/bin/sh

#Script to download records from nucleotide database

NUCLEOTIDE_DATABASE="nucleotide"
QUERY="(human[organism] OR \"homo sapiens\"[organism]) AND (\"mitochondrial\"[title] or mitochondrion[TITLE]) AND \"complete genome\"[TITLE] "

RESULT=$(esearch -db "$NUCLEOTIDE_DATABASE" -query "$QUERY" | epost -format xml)
echo "$RESULT"
WEBENV=$(echo "$RESULT" | xtract -pattern ENTREZ_DIRECT -element WebEnv)
QUERY_KEY=$(echo "$RESULT" | xtract -pattern ENTREZ_DIRECT -element QueryKey)
COUNT=$(echo "$RESULT" | xtract -pattern ENTREZ_DIRECT -element Count)

echo "Total Records Found: $COUNT"
echo "WebEnv: $WEBENV"

# Loop through results in chunks of 500
CHUNK_SIZE=500

for i in $(seq 0 $CHUNK_SIZE $COUNT); do
    echo "Fetching records $i to $((i + CHUNK_SIZE))..."

    curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$NUCLEOTIDE_DATABASE&query_key=$QUERY_KEY&WebEnv=$WEBENV&retstart=$i&retmax=$CHUNK_SIZE&rettype=gb&retmode=xml" \
        >>nucleotide_records_test.xml

done

echo "All records fetched successfully!"
