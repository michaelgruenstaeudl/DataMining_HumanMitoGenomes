# NCBI Record Mining

Files and scripts to identify associated entries across NCBI GenBank, NCBI SRA, and other NCBI databases

## FILE CONTENTS

### Code and Scripts

- `CODE_ncbi_nucleotide_fetcher.sh`: A Bash script to search for and then download records from NCBI Nucleotide in XML format.

- `CODE_Comparison_Mitobank_vs_EntrezQuery.ipynb`: An IPython Notebook that lists various computational comparisons between the GenBank records listed on Mitobank and those extracted via Entrez Direct, including comparisons of the number of records and which ones are unique to each file.

- `CODE_EntrezQueries_NumberOfRecordsInDifferentNCBIDatabases.md`: A markdown file that lists various Entrez Direct queries and their results.

- `CODE_nucleotide_sra_data_extraction.ipynb`: An IPython Notebook that contains scripts for extracting nucleotide metadata in tabular format.

### Data Files

- `DATA_GenBank_accessions_Mitobank_2025_01_31.txt`: A text file listing the GenBank accession numbers of all full-length human mitochondrial genome sequences listed on Mitobank as of 31-Jan-2024; it was downloaded directly from [Mitobank](https://www.mitomap.org/foswiki/bin/view/MITOMAP/Mitobank). The file contains a total of 61845 accession numbers.

- `DATA_GenBank_accessions_EntrezQuery_2025_01_31.txt`: A text file listing the GenBank accession numbers of all complete human mitochondrial genome sequences received when running an Entrez Direct query as of 31-Jan-2024. Section _Downloaded accession IDs from nucleotide database_ of file `entrez_scripts_and_results.md` specifies the Entrez Direct query used. The file contains a total of 62167 accession numbers.

- `DATA_Nucleotide_Metadata.csv`: A comma-separated list of all nucleotide metadata records extracted using script `extract_nucleotide_metadata.ipynb`. Accession Id, BioProject Id, BioSample Id, Publication Title, Publication Reference.

- `DATA_Nucleotide_Summary_records.csv`: A comma-separated list of all nucleotide summary records extracted using the following query: `"Homo sapiens[ORGN] AND complete genome[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] NOT (unverified OR Homo sp. Altai OR Denisova hominin OR neanderthalensis OR heidelbergensis OR consensus)"`

## RESOURCES

### Obtaining publications of NCBI Nucleotide records
Bash code to extract the publications that the 62,173 complete human mitogenomes stored on NCBI Nucleotide were published in:
`cat DATA_Nucleotide_Metadata.csv | xsv select 9-11 | grep -v "Unpublished\|Submitted" | sort -u`


### Promising Entrez search strategies
Up-to-date list of **all GenBank accession numbers** of all **complete mitochondrial genomes of humans** (total: 61,845) stored on GenBank: https://www.mitomap.org/cgi-bin/genbank_ids.cgi

#### On NCBI GenBank

`Homo sapiens[ORGN] AND complete genome[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] NOT (unverified OR Homo sp. Altai OR Denisova hominin OR neanderthalensis OR heidelbergensis OR consensus)` produces 62,173 hits on NCBI Nucleotide

#### On NCBI SRA

`(mitochondrion[ALL] OR mitochondrial[ALL]) AND ("Homo sapiens"[ORGN] OR human[TITLE]) AND "biomol dna"[PROP] AND "platform illumina"[PROP]` produces 14,238 hits on NCBI SRA
