#!/usr/bin/env nextflow

process mitoz_assembler {

    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_file_list), val(is_single_end)

    output:
    tuple val(sample_id), path("megahit/mito_denovo.megahit.result/mito_denovo.megahit.mitogenome.fa"), emit: mitoz_assembler_output
    tuple val(sample_id), path("megahit"), emit: mitoz_assembler_publish_dir

    script:
    def input_file_param = is_single_end ? "--fq1 ${input_file_list[0]}" : "--fq1 ${input_file_list[0]} --fq2 ${input_file_list[1]}"
    """
    #Env variable for mitoz
    export PATH=/app/anaconda/envs/mitoz3.6/bin:\$PATH

    # Run MitoZ assembler
    mitoz assemble   \
        --genetic_code 2   \
        --clade Chordata   \
        --outprefix mito_denovo   \
        --thread_number 16   \
        ${input_file_param} \
        --assembler megahit   \
        --kmers_megahit 21 29 39 59 79 99 119 141 \
        --requiring_taxa 9606
    """
}
