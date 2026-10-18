#!/bin/bash

#Script to count and view mtDNA nucleotide records in NCBI Nucleotide database using Entrez Direct (EDirect) tools

NUCLEOTIDE_DATABASE="nucleotide"
QUERY="(human[organism] OR \"homo sapiens\"[organism]) AND (\"mitochondrial\"[title] or mitochondrion[TITLE]) AND \"complete genome\"[TITLE] "

esearch -db "$NUCLEOTIDE_DATABASE" -query "$QUERY"

# Script to count and view human mtDNA SRA records available in SRA database
esearch \
    -db SRA \
    -query "(human[organism] OR \"homo sapiens\"[organism]) \
    AND (\"mitochondrial\"[title] or mitochondrion[TITLE])"
