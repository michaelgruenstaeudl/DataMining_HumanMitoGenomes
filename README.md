# NCBI Record Mining
Scripts to identify associated entries across NCBI GenBank, NCBI SRA, and other NCBI databases

## FILE CONTENTS
- `2000until2010_GenBank_accessions_plus_metadata.tsv`: a tab-separated list of all GenBank record accessions between 2000 and 2010 with valuable metadata extracted. For instance, with `awk -F'\t' '{print$9}' 2000until2010_GenBank_accessions_plus_metadata.tsv | sort -u | grep -v "Submitted"` you can extract all the publications that describe these GenBank records; in these publications, the SRA numbers may be saved.

## RESOURCES
Up-to-date list of **all GenBank accession numbers** of all **complete mitochondrial genomes of humans** (total: 61,845) stored on GenBank: https://www.mitomap.org/cgi-bin/genbank_ids.cgi

### Promising Entrez search strategies

#### On NCBI GenBank
```Homo sapiens[ORGN] AND complete genome[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] NOT (unverified OR Homo sp. Altai OR Denisova hominin OR neanderthalensis OR heidelbergensis OR consensus)``` produces 62,173 hits on NCBI Nucleotide
#### On NCBI SRA
```(mitochondrion[ALL] OR mitochondrial[ALL]) AND ("Homo sapiens"[ORGN] OR human[TITLE]) AND "biomol dna"[PROP] AND "platform illumina"[PROP]``` produces 14,238 hits on NCBI SRA


