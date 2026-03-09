### How to circularize a human mitogenome
Human mitochondrial genome assemblies often appear linear with duplicated ends because assemblers do not enforce a circular topology.
The following steps (a) rotate the sequence so all genomes start at the D-loop (control region) and (b) trim duplicated 3' region or expand an incomplete 3' region.

The input files are located in folder 'circularlizing_mitogenomes'

#### 1. Using rotate to rotate mitogenome to correct start
https://github.com/richarddurbin/rotate
See examples in: https://wellcomeopenresearch.org/articles/8-401

##### Installation
```bash
$ git clone https://github.com/richarddurbin/rotate.git ; cd rotate ; make
```

##### Usage: rotate to anchor string based on first 50 bp of D-loop, allowing for 5 mismatches

###### Mitogenome has wrong start and is too long
```bash
$ ./rotate -s GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCAT -m 5 MG936619_too_long_and_wrong_start.fasta > MG936619_too_long.fasta

# If more than two matches detected, take the latter match position
sequence        MG936619_too_long_and_wrong_start       length  18465
  forward       8299    nerr    0
  forward       10196   nerr    0

$ ./rotate -x 10196 MG936619_too_long_and_wrong_start.fasta > MG936619_too_long.fasta
$ sed -i '1s/.*/>MG936619_too_long__startcorrected/' MG936619_too_long.fasta
```
###### Mitogenome has wrong start and is too short
```bash
$ ./rotate -s GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCAT -m 5 MG936619_too_short_and_wrong_start.fasta > MG936619_too_short.fasta
$ sed -i '1s/.*/>MG936619_too_short__startcorrected/' MG936619_too_short.fasta
```

#### 2a. If rotated mitogenome too long, trim off extra end
The following Python script `trim_terminal_suffix_dup.py` reads a FASTA sequence, finds the longest prefix of the sequence that is also an identical suffix, and removes that duplicated suffix from the 3' end, printing the trimmed sequence only if such an overlap exists.

```python
#!/usr/bin/env python3

import sys

if len(sys.argv) != 2:
    sys.exit("Usage: trim_terminal_prefix_dup.py input.fasta")

fasta = sys.argv[1]

header = None
seq_parts = []

with open(fasta) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            header = line
        else:
            seq_parts.append(line)

seq = "".join(seq_parts)
n = len(seq)

# Compute KMP prefix function
pi = [0] * n
j = 0

for i in range(1, n):
    while j > 0 and seq[i] != seq[j]:
        j = pi[j-1]
    if seq[i] == seq[j]:
        j += 1
    pi[i] = j

overlap = pi[-1]

# Print only if trimming occurs
if overlap > 0:
    trimmed = seq[:-overlap]
    new_header = header + "_suffix_dup_trimmed"
    print(new_header)
    print(trimmed)
else:
    sys.exit("No overlap")
```

```bash
python trim_terminal_suffix_dup.py MG936619_too_long.fasta > MG936619_suffix_dup_trimmed.fasta
```


#### 2b. If rotated mitogenome too short, extend via Novoplasty
Use "partial_mito.fasta" as "seed input":
```
Project:
-----------------------
Project name          = mito_extend
Type                  = mito
Genome Range          = 14000-20000
K-mer                 = 39
Max memory            =
Extended log          = 0
Save assembled reads  = yes
Seed Input            = partial_mito.fasta
Reference sequence    =
Variance detection    =
Chloroplast sequence  =

Dataset 1:
-----------------------
Read Length           = 150
Insert size           = 300
Platform              = illumina
Single/Paired         = PE
Combined reads        =
Forward reads         = reads_R1.fastq.gz
Reverse reads         = reads_R2.fastq.gz

Optional:
-----------------------
Insert size auto      = yes
Use Quality Scores    = yes
```

#### If 2b. successful, do 2a.
