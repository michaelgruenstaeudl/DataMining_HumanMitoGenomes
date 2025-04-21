#!/usr/bin/env nextflow

// container 'genomicpariscentre/trimgalore:0.6.10' 

process qualityControl {
    publishDir 'results', mode: 'copy'
    container 'community.wave.seqera.io/library/trim-galore:0.6.10--1bf8ca4e1967cd18'

    input:
    path input_file1
    path input_file2

    output:
    path "trimmed_output/*"

    script:
    """
    echo Input file path: ${input_file1} 
    trim_galore --paired --length 30 --quality 20 \
    --fastqc \
    --output_dir trimmed_output \
    ${input_file1} ${input_file2}
    """
}
