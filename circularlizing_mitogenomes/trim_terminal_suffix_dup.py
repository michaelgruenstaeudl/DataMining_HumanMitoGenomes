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
    new_header = header + "_suffix_dupl_trimmed"
    print(new_header)
    print(trimmed)
else:
    sys.exit("No overlap")
