# NCBI Record Mining
Scripts to identify associated entries across NCBI GenBank, NCBI SRA, and other NCBI databases

### Resources
Up-to-date list of **all GenBank accession numbers** of all **complete mitochondrial genomes of humans** (total: 61,845) stored on GenBank: https://www.mitomap.org/cgi-bin/genbank_ids.cgi

### Promising Entrez search strategies

#### On NCBI GenBank
```Homo sapiens[ORGN] AND complete genome[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] NOT (unverified OR Homo sp. Altai OR Denisova hominin OR neanderthalensis OR heidelbergensis OR consensus)``` produces 62,173 hits on NCBI Nucleotide
#### On NCBI SRA
```(mitochondrion[ALL] OR mitochondrial[ALL]) AND ("Homo sapiens"[ORGN] OR human[TITLE]) AND "biomol dna"[PROP] AND "platform illumina"[PROP]``` produces 14,238 hits on NCBI SRA
