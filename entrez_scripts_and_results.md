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
