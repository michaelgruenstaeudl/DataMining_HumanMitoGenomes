# NCBI Record Mining

Files and scripts to identify associated entries across NCBI GenBank, NCBI SRA, and other NCBI databases

## FILE CONTENTS

- **2000until2010_GenBank_accessions_plus_metadata.tsv**:
  A tab-separated list of all GenBank record accessions between 2000 and 2010 with valuable metadata extracted. For instance, with `awk -F'\t' '{print$9}' 2000until2010_GenBank_accessions_plus_metadata.tsv | sort -u | grep -v "Submitted"` you can extract all the publications that describe these GenBank records; in these publications, the SRA numbers may be saved.

- **ncbi_nucleotide_fetcher.sh**
  A Bash script to download GenBank records from Nucleotide database in XML format.

- **genbank_ids.txt**:
  A text file listing the GenBank accession numbers of all full-length human mitochondrial genome sequences listed on Mitobank as of 31-Jan-2024; it was downloaded directly from [Mitobank](https://www.mitomap.org/foswiki/bin/view/MITOMAP/Mitobank). The file contains a total of 61845 accession numbers.

- **homosapian_nucleotide_accession_list.txt**:
  A text file listing the GenBank accession numbers of all complete human mitochondrial genome sequences received when running an Entrez Direct query as of 31-Jan-2024. Section *Downloaded accession IDs from nucleotide database* of file `entrez_scripts_and_results.md` specifies the Entrez Direct query used. The file contains a total of 62167 accession numbers.

- **mitochondrial_genome_analysis.ipynb**:
  An IPython Notebook that lists various computational analyses of the GenBank records listed in `genbank_ids.txt` and `homosapian_nucleotide_accession_list.txt` (e.g., comparing the common records between them and identifying the unique records of each file).

## RESOURCES

Up-to-date list of **all GenBank accession numbers** of all **complete mitochondrial genomes of humans** (total: 61,845) stored on GenBank: https://www.mitomap.org/cgi-bin/genbank_ids.cgi

### Promising Entrez search strategies

#### On NCBI GenBank

`Homo sapiens[ORGN] AND complete genome[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] NOT (unverified OR Homo sp. Altai OR Denisova hominin OR neanderthalensis OR heidelbergensis OR consensus)` produces 62,173 hits on NCBI Nucleotide

#### On NCBI SRA

`(mitochondrion[ALL] OR mitochondrial[ALL]) AND ("Homo sapiens"[ORGN] OR human[TITLE]) AND "biomol dna"[PROP] AND "platform illumina"[PROP]` produces 14,238 hits on NCBI SRA
