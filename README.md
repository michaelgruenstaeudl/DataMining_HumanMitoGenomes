# NCBI Record Mining
Scripts to identify associated entries across NCBI GenBank, NCBI SRA, and other NCBI databases

### Resources
Up-to-date list of *all GenBank accession numbers* of all *complete mitochondrial genomes of humans* stored on GenBank:
https://www.mitomap.org/cgi-bin/genbank_ids.cgi

### Promising search strategies

#### On SRA
```(mitochondrion[ALL] OR mitochondrial[ALL]) AND ("Homo sapiens"[ORGN] OR human[TITLE]) AND "biomol dna"[PROP] AND "platform illumina"[PROP]``` produces 14238 hits
