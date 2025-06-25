# Data Mining on Human Mitochondrial Genome Data

Pipeline for data mining on human mitochondrial genome data

## FILE CONTENTS

### Code and Scripts

- `CODE_ncbi_nucleotide_fetcher.sh`: A Bash script to search for and then download records from NCBI Nucleotide in XML format.

- `CODE_Comparison_Mitobank_vs_EntrezQuery.ipynb`: An IPython Notebook that lists various computational comparisons between the GenBank records listed on Mitobank and those extracted via Entrez Direct, including comparisons of the number of records and which ones are unique to each file.

- `CODE_EntrezQueries_NumberOfRecordsInDifferentNCBIDatabases.md`: A markdown file that lists various Entrez Direct queries and their results.

- `CODE_nucleotide_sra_data_extraction.ipynb`: An IPython Notebook that contains scripts for extracting nucleotide metadata in tabular format.

- `CODE_SRA_metadata_extractor.sh`: A shell script to extract SRA metadata records related to human mitochondrial nucleotide records through biosample and bioproject database.

- `CODE_Nucleotide_SRA_data_analysis.ipynb`: An IPython Notebook that contains analysis of nucleotide metadata with SRA metadata based on BioProject Id and BioSampleId.

- `CODE_mitochondrial_pubmed_article_extraction.py`: A Python Script to extract pubmed record informations from the list of title in CSV file with header `TITLE'.

  ```
  For example:
  |   TITLE  |
  |----------|
  | title_a  |
  | title_b  |
  | title_c  |
  ```

  Command-line arguments containing:

  - mail (str): Email address for PubMed API.
  - verbose (bool): Verbosity flag for logging.
  - filepath (str): Path to the input CSV file containing article titles.

  `-m "example@email.com" -f "DATA_publications_of_human_mitogenomes_on_NCBINulceotide.csv"`

  This script generates following documents:

  - Extract and save PubMed metadata for each title in the CSV file in `DATA_pubmed_metadata.csv`.
  - Extract full text from PubMed articles and store in `DATA_pubmed_records.json`.
  - Save the pubmed records which contains data source informations in `DATA_pubmed_records_with_data_source_info.json`.

- `log` directory: This directory contains log generated from `CODE_mitochondrial_pubmed_article_extraction.py` script execution.

- `CODE_mitochondrial_sra_info_extraction.py`: A python script to extract the SRA record informations and first paragraph/data from corresponding pubmed article based on pubmed article title, BioProject Accession Number or European Nucleotide Archive Number or SRA ID or BioSample Accession Number.

### Data Files

- `DATA_GenBank_accessions_Mitobank_2025_01_31.txt`: A text file listing the GenBank accession numbers of all full-length human mitochondrial genome sequences listed on Mitobank as of 31-Jan-2024; it was downloaded directly from [Mitobank](https://www.mitomap.org/foswiki/bin/view/MITOMAP/Mitobank). The file contains a total of 61845 accession numbers.

- `DATA_GenBank_accessions_EntrezQuery_2025_01_31.txt`: A text file listing the GenBank accession numbers of all complete human mitochondrial genome sequences received when running an Entrez Direct query as of 31-Jan-2024. Section _Downloaded accession IDs from nucleotide database_ of file `entrez_scripts_and_results.md` specifies the Entrez Direct query used. The file contains a total of 62167 accession numbers.

- `DATA_Nucleotide_Metadata.csv`: A comma-separated list of all nucleotide metadata records extracted using script `extract_nucleotide_metadata.ipynb`. Accession Id, BioProject Id, BioSample Id, Publication Title, Publication Reference.

- `DATA_Nucleotide_Summary_records.csv`: A comma-separated list of all nucleotide summary records extracted using the following query: `"Homo sapiens[ORGN] AND complete genome[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] NOT (unverified OR Homo sp. Altai OR Denisova hominin OR neanderthalensis OR heidelbergensis OR consensus)"`

- `DATA_SRA_Metadata_linked_through_BioSample.csv`: A comma-separated list of all SRA metadata records linked to human mitochondrial nucleotide records through BioSample database.

- `DATA_SRA_Metadata_linked_through_BioProject.csv`:A comma-separated list of all SRA metadata records linked to human mitochondrial nucleotide records through BioProject database.

- `DATA_pubmed_metadata.csv`: A comma-separated list of all pubmed metadata records linked to human mitochondrial nucleotide records. This CSV contains following informations: "Pubmed_ID", "Title", "Published_Year", "DataBankList", "Full_Article_URL", "is_PMC", "Error".

- `DATA_pubmed_records_with_data_source_info.json`: A JSON file contains all pubmed article informations which contains data source informations.

- `pubmed_records.json`: A JSON file containing pubmed article informations.
- `DATA_pubmed_info_with_bioproj_and_sra.csv`: A CSV file containing pubmed article with corresponding BioProject Id or European Nucleotide Archive ID or SRA ID. This file is important for `CODE_mitochondrial_sra_info_extraction.py`.
- `DATA_pubmed_info_with_bioproj_sra_sample.csv`: A CSV file containing following informations:

  - Title: Title of Pubmed Article
  - BioProject_uid: BioProject Unique Identifier associated to the pubmed article
  - BioProject: BioProject Accession Number associated to the pubmed article
  - SRA_Code : SRA code associated with pubmed Article
  - Pubmed_ID : Pubmed Unique Identifier related to article.
  - PMC_ID : Pubmed central unique identifier (Pubmed article id)
  - SRA_Count: Total number of SRA records associated with BioProject.
  - SRA_Id_list: Unique SRA Id list associated with BioProject.
  - Is_in_Genbank: Whether records present in NCBI Genbank or European Nucleotide Archive.

- `DATA_sra_info_in_pmc.json`: A JSON file which contains SRA records metadata. It contains all SRA records present in `DATA_pubmed_info_with_bioproj_sra_sample.csv`.

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

## Funding & Citation
This code development is part of NIH project R01LM014506.
