#!/usr/bin/env nextflow

process bbmap_index {

    label "bbmap_env"
    publishDir "${parent_output_dir}", mode: 'copy', overwrite: true

    input:
    path reference_fasta
    val parent_output_dir

    output:

    path "ref_index"

    script:
    """
    mkdir ref_index
    bbmap.sh \
        ref=${reference_fasta} \
        path=ref_index
    """
}
