#!/usr/bin/env nextflow


process qualityControl {
    publishDir 'results', mode: 'copy'

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
