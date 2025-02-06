# NCBI Record Mining

Files and scripts to identify associated entries across NCBI GenBank, NCBI SRA, and other NCBI databases

## FILE CONTENTS

### Scripts

- **ncbi_nucleotide_fetcher.sh**: A Bash script to search for and then download records from NCBI Nucleotide in XML format.

- **Comparison_Mitobank_vs_EntrezQuery.ipynb**: An IPython Notebook that lists various computational analyses of the GenBank records listed on Mitobank versus those extracted via Entrez Direct. The analyses include, for instance, a comparison of the number of records and which ones are unique to each file.

- **EntrezQueries_and_results**STEP01**NumberOfRecordsInDifferentNCBIDatabases.md**: A markdown file that lists various Entrez Direct queries and their results.

- **extract_nucleotide_metadata.ipynb**: An IPython Notebook that contains scripts for extracting nucleotide metadata in tabuler format.

### Data Files

- **GenBank_accessions**Mitobank**2025_01_31.txt**: A text file listing the GenBank accession numbers of all full-length human mitochondrial genome sequences listed on Mitobank as of 31-Jan-2024; it was downloaded directly from [Mitobank](https://www.mitomap.org/foswiki/bin/view/MITOMAP/Mitobank). The file contains a total of 61845 accession numbers.

- **GenBank_accessions**EntrezQuery**2025_01_31.txt**: A text file listing the GenBank accession numbers of all complete human mitochondrial genome sequences received when running an Entrez Direct query as of 31-Jan-2024. Section _Downloaded accession IDs from nucleotide database_ of file `entrez_scripts_and_results.md` specifies the Entrez Direct query used. The file contains a total of 62167 accession numbers.

- **2000until2010_GenBank_accessions_plus_metadata.tsv**: A tab-separated list of all GenBank record submitted between 2000 and 2010 including valuable metadata. For instance, the Bash command `awk -F'\t' '{print$9}' 2000until2010_GenBank_accessions_plus_metadata.tsv | sort -u | grep -v "Submitted"` extracts all the publications that describe these GenBank records; in these publications, the SRA numbers may be located.

- **Nucleotide_Metadata.csv**: A comma-separated list of all nucleotide metadata records extracted using **extract_nucleotide_metadata.ipynb** script. Accession Id, BioProject Id, BioSample Id, Title, References.

- **Nucleotide_Summary_records.csv**: A comma-separated list of all nucleotide summary records extracted using following query: `"Homo sapiens[ORGN] AND complete genome[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] NOT (unverified OR Homo sp. Altai OR Denisova hominin OR neanderthalensis OR heidelbergensis OR consensus)"`

## RESOURCES

Up-to-date list of **all GenBank accession numbers** of all **complete mitochondrial genomes of humans** (total: 61,845) stored on GenBank: https://www.mitomap.org/cgi-bin/genbank_ids.cgi

### Promising Entrez search strategies

#### On NCBI GenBank

`Homo sapiens[ORGN] AND complete genome[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] NOT (unverified OR Homo sp. Altai OR Denisova hominin OR neanderthalensis OR heidelbergensis OR consensus)` produces 62,173 hits on NCBI Nucleotide

#### On NCBI SRA

`(mitochondrion[ALL] OR mitochondrial[ALL]) AND ("Homo sapiens"[ORGN] OR human[TITLE]) AND "biomol dna"[PROP] AND "platform illumina"[PROP]` produces 14,238 hits on NCBI SRA
