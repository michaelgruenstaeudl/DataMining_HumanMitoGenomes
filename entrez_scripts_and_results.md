# Entrez Direct Scripts and Results

This document contains shell scripts and their outputs for analysis purposes, including the date and time of execution. It serves as a comprehensive history of the scripts executed and tracks their outputs.

## Total Number of Complete Human Mitochondrial Genomes in Nucleotide Database

**Execution Date**: January 28, 2025

**Script**:

```
esearch \
    -db nucleotide \
    -query "(human[organism] OR \"homo sapiens\"[organism]) \
    AND (\"mitochondrial\"[title] or mitochondrion[TITLE]) \
    AND complete genome[TITLE] "
```

**Result**:

```
<ENTREZ_DIRECT>
  <Db>nuccore</Db>
  <WebEnv>MCID_67a268ba532c973c850e2703</WebEnv>
  <QueryKey>1</QueryKey>
  <Count>62753</Count>
  <Step>1</Step>
  <Elapsed>1</Elapsed>
</ENTREZ_DIRECT>
```

## Total Number of records in BioProject Database related to Complete Human Mitochondrial Genomes

**Execution Date**: January 28, 2025

**Script**:

```
esearch \
    -db bioproject \
    -query "(human[organism] OR \"homo sapiens\"[organism]) \
    AND (\"mitochondrial\"[title] or mitochondrion[TITLE]) \
    AND complete genome[TITLE] "
```

**Result**:

```
<ENTREZ_DIRECT>
  <Db>bioproject</Db>
  <WebEnv>MCID_67a26dcada8b40275c0db956</WebEnv>
  <QueryKey>1</QueryKey>
  <Count>0</Count>
  <Step>1</Step>
  <Elapsed>1</Elapsed>
</ENTREZ_DIRECT>
```

## Total Number of records in BioSample Database related to Complete Human Mitochondrial Genomes

**Execution Date**: January 28, 2025

**Script**:

```
esearch \
    -db BioSample \
    -query "(human[organism] OR \"homo sapiens\"[organism]) \
    AND (\"mitochondrial\"[title] or mitochondrion[TITLE]) \
    AND complete genome[TITLE] "
```

**Result**:

```
<ENTREZ_DIRECT>
  <Db>biosample</Db>
  <WebEnv>MCID_67a26e0209a1154b8c090595</WebEnv>
  <QueryKey>1</QueryKey>
  <Count>452</Count>
  <Step>1</Step>
  <Elapsed>1</Elapsed>
</ENTREZ_DIRECT>
```

## Total Number of records in SRA Database related to Human Mitochondrial Genomes

**Execution Date**: January 29, 2025

**Script**:

```
esearch \
    -db SRA \
    -query "(human[organism] OR \"homo sapiens\"[organism]) \
    AND (\"mitochondrial\"[title] or mitochondrion[TITLE])"
```

**Result**:

```
<ENTREZ_DIRECT>
  <Db>sra</Db>
  <WebEnv>MCID_67a26e67db4b35cd8c0683af</WebEnv>
  <QueryKey>1</QueryKey>
  <Count>25831</Count>
  <Step>1</Step>
  <Elapsed>1</Elapsed>
</ENTREZ_DIRECT>
```

## Total Number of records in SRA Database related to Human Mitochondrial Genomes and having biomol DNA and platform illumina properties.

**Execution Date**: January 29, 2025

**Script**:

```
esearch \
    -db SRA \
    -query "(\"Homo sapiens\"[ORGN] OR human[TITLE]) \
    AND (mitochondrion[ALL] OR mitochondrial[ALL]) \
    AND \"biomol dna\"[PROP] AND \"platform illumina\"[PROP]"
```

**Result**:

```
<ENTREZ_DIRECT>
  <Db>sra</Db>
  <WebEnv>MCID_67a274e6e1c24557800a865b</WebEnv>
  <QueryKey>1</QueryKey>
  <Count>14361</Count>
  <Step>1</Step>
  <Elapsed>1</Elapsed>
</ENTREZ_DIRECT>
```

## Total Number of records in Nucleotide Database related to Human Mitochondrial Genomes and having sequence length range from 15400 to 16700

**Execution Date**: January 29, 2025

**Script**:

```
esearch \
    -db nucleotide \
    -query "\"Homo sapiens\"[ORGN] AND \"complete genome\"[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] \
    NOT (unverified OR \"Homo sp. Altai\" OR \"Denisova hominin\" OR neanderthalensis OR heidelbergensis OR consensus)"
```

**Result**:

```
<ENTREZ_DIRECT>
  <Db>nuccore</Db>
  <WebEnv>MCID_67a276bb2fb23d4667016b90</WebEnv>
  <QueryKey>1</QueryKey>
  <Count>62173</Count>
  <Step>1</Step>
  <Elapsed>1</Elapsed>
</ENTREZ_DIRECT>
```

## Total Number of records in BioProject Database linked with Nucleotide database related to Complete Human Mitochondrial Genomes

**Execution Date**: January 30, 2025

**Script**:

```
esearch \
    -db nucleotide \
    -query "\"Homo sapiens\"[ORGN] AND \"complete genome\"[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] \
    NOT (unverified OR \"Homo sp. Altai\" OR \"Denisova hominin\" OR neanderthalensis OR heidelbergensis OR consensus)" |
    elink -target bioproject
```

**Result**:

```
<ENTREZ_DIRECT>
  <Db>bioproject</Db>
  <WebEnv>MCID_67a289d20959ac373907af17</WebEnv>
  <QueryKey>102</QueryKey>
  <Count>5</Count>
  <Step>2</Step>
  <Elapsed>125</Elapsed>
</ENTREZ_DIRECT>
```

## Total Number of records in BioSample Database linked with Nucleotide database related to Complete Human Mitochondrial Genomes

**Execution Date**: January 30, 2025

**Script**:

```
esearch \
    -db nucleotide \
    -query "\"Homo sapiens\"[ORGN] AND \"complete genome\"[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] \
    NOT (unverified OR \"Homo sp. Altai\" OR \"Denisova hominin\" OR neanderthalensis OR heidelbergensis OR consensus)" |
    elink -target BioSample
```

**Result**:

```
<ENTREZ_DIRECT>
  <Db>BioSample</Db>
  <WebEnv>MCID_67a282004a648e9a3306a676</WebEnv>
  <QueryKey>102</QueryKey>
  <Count>15</Count>
  <Step>2</Step>
  <Elapsed>87</Elapsed>
</ENTREZ_DIRECT>
```

## Total Number of records in SRA Database linked with Nucleotide database related to Complete Human Mitochondrial Genomes

**Execution Date**: January 31, 2025

**Script**:

```
esearch \
    -db nucleotide \
    -query "\"Homo sapiens\"[ORGN] AND \"complete genome\"[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] \
    NOT (unverified OR \"Homo sp. Altai\" OR \"Denisova hominin\" OR neanderthalensis OR heidelbergensis OR consensus)" |
    elink -target SRA
```

**Result**:

```
<ENTREZ_DIRECT>
  <Db>sra</Db>
  <WebEnv>MCID_67a28ab42db8ab41440d0100</WebEnv>
  <QueryKey>1</QueryKey>
  <Count>0</Count>
  <Step>2</Step>
  <Elapsed>482</Elapsed>
</ENTREZ_DIRECT>
```

## Total Number of records in SRA Database linked with Nucleotide database through BioProject database related to Complete Human Mitochondrial Genomes (Nucleotide -> BioProject -> SRA)

**Execution Date**: January 31, 2025

**Script**:

```
esearch \
    -db nucleotide \
    -query "\"Homo sapiens\"[ORGN] AND \"complete genome\"[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] \
    NOT (unverified OR \"Homo sp. Altai\" OR \"Denisova hominin\" OR neanderthalensis OR heidelbergensis OR consensus)" |
    elink -target bioproject |
    elink -target SRA
```

**Result**:

```
<ENTREZ_DIRECT>
  <Db>SRA</Db>
  <WebEnv>MCID_67a2bfa78d5f644c7503631e</WebEnv>
  <QueryKey>104</QueryKey>
  <Count>259</Count>
  <Step>3</Step>
  <Elapsed>122</Elapsed>
</ENTREZ_DIRECT>
```

With different search condition: \

**Script**:

```
esearch -db nucleotide -query "(human[organism] OR \"homo sapiens\"[organism]) \
    AND (\"mitochondrial\"[title] or mitochondrion[TITLE]) \
    AND complete genome[TITLE]" |
    elink -target bioproject |
    elink -target SRA
```

**Result**:

```
<ENTREZ_DIRECT>
  <Db>SRA</Db>
  <WebEnv>MCID_67a2dbe000592acfb406c179</WebEnv>
  <QueryKey>104</QueryKey>
  <Count>277</Count>
  <Step>3</Step>
  <Elapsed>261</Elapsed>
</ENTREZ_DIRECT>
```

## Total Number of records in SRA Database linked with Nucleotide database through BioSample database related to Complete Human Mitochondrial Genomes (Nucleotide -> BioSample -> SRA)

**Execution Date**: January 31, 2025

**Script**:

```
esearch \
    -db nucleotide \
    -query "\"Homo sapiens\"[ORGN] AND \"complete genome\"[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] \
    NOT (unverified OR \"Homo sp. Altai\" OR \"Denisova hominin\" OR neanderthalensis OR heidelbergensis OR consensus)" |
    elink -target biosample |
    elink -target SRA
```

**Result**:

```
<ENTREZ_DIRECT>
  <Db>SRA</Db>
  <WebEnv>MCID_67a2c0a17ec82cf03b058e81</WebEnv>
  <QueryKey>104</QueryKey>
  <Count>744</Count>
  <Step>3</Step>
  <Elapsed>127</Elapsed>
</ENTREZ_DIRECT>
```

With different search condition: \

**Script**:

```
esearch -db nucleotide -query "(human[organism] OR \"homo sapiens\"[organism]) \
    AND (\"mitochondrial\"[title] or mitochondrion[TITLE]) \
    AND complete genome[TITLE]" |
    elink -target biosample |
    elink -target SRA
```

**Result**:

```
<ENTREZ_DIRECT>
  <Db>SRA</Db>
  <WebEnv>MCID_67a2c60ffb3ac7f23d0f7b85</WebEnv>
  <QueryKey>104</QueryKey>
  <Count>746</Count>
  <Step>3</Step>
  <Elapsed>260</Elapsed>
</ENTREZ_DIRECT>
```

## Total Number of records in SRA Database linked with Nucleotide database through BioSample and BioProject database related to Complete Human Mitochondrial Genomes (Nucleotide -> BioSample -> BioProject -> SRA)

**Execution Date**: February 3, 2025

**Script**:

```
esearch \
    -db nucleotide \
    -query "\"Homo sapiens\"[ORGN] AND \"complete genome\"[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] \
    NOT (unverified OR \"Homo sp. Altai\" OR \"Denisova hominin\" OR neanderthalensis OR heidelbergensis OR consensus)" |
    elink -target biosample |
    elink -target bioproject |
    elink -target SRA
```

**Result**:

```
<ENTREZ_DIRECT>
  <Db>SRA</Db>
  <WebEnv>MCID_67a2e1b12a74d727b60da31f</WebEnv>
  <QueryKey>106</QueryKey>
  <Count>2843</Count>
  <Step>4</Step>
  <Elapsed>178</Elapsed>
</ENTREZ_DIRECT>
```

With another query

**Script**:

```
esearch \
    -db nucleotide \
    -query "(human[organism] OR \"homo sapiens\"[organism]) \
    AND (\"mitochondrial\"[title] or mitochondrion[TITLE]) \
    AND complete genome[TITLE] " |
    elink -target BioSample |
    elink -target BioProject |
    elink -target SRA
```

**Result**:

```
<ENTREZ_DIRECT>
  <Db>SRA</Db>
  <WebEnv>MCID_67a2e36ed12ee4095d0e634b</WebEnv>
  <QueryKey>106</QueryKey>
  <Count>1119</Count>
  <Step>4</Step>
  <Elapsed>89</Elapsed>
</ENTREZ_DIRECT>
```

## Downloaded accession ids from nucleotide database

**Execution Date**: February 4, 2025

**Script**:

```
esearch -db nucleotide -query "\"Homo sapiens\"[ORGN] \
    AND \"complete genome\"[TITLE] AND mitochondrion[FILT] AND 015400:016700[SLEN] \
    NOT (unverified OR \"Homo sp. Altai\" OR \"Denisova hominin\" OR neanderthalensis OR heidelbergensis OR consensus)" |
    efetch -format acc -mode text >homosapian_nucleotide_accession_list.txt
```

**Result**:

- 62173 accession ids are downloaded and stored in text file _homosapian_nucleotide_accession_list.txt_.
- We also downloaded full length mitochondrial accession id lists from mitobank. Filename: _genbank_ids.txt_
- Link to mitobank: https://www.mitomap.org/foswiki/bin/view/MITOMAP/Mitobank

```
<ENTREZ_DIRECT>
  <Db>nuccore</Db>
  <WebEnv>MCID_67a41b51643a9255a30f8ad7</WebEnv>
  <QueryKey>1</QueryKey>
  <Count>62173</Count>
  <Step>1</Step>
  <Elapsed>1</Elapsed>
</ENTREZ_DIRECT>
```
