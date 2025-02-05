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

## Total Number of records in BioProject Database linked with Nucleotide database related Complete Human Mitochondrial Genomes

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

## Total Number of records in BioSample Database linked with Nucleotide database related Complete Human Mitochondrial Genomes

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

## Total Number of records in SRA Database linked with Nucleotide database related Complete Human Mitochondrial Genomes

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
