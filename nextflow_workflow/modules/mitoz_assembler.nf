#!/usr/bin/env nextflow

// output:
//     path "novoplasty_outptut/*"
process mitoz_assembler {
    // container 'apptainer/MitoZ_v3.6.sif'

    tag "${sample_id}"

    publishDir "${parent_output_dir}/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(input_file1), path(input_file2)
    val parent_output_dir
    val projectDir
    val eteDbDir

    output:
    path "megahit/*"

    script:
    """
    #Env variable for mitoz
    export PATH=/app/anaconda/envs/mitoz3.6/bin:\$PATH

    # Run MitoZ assembler
    mitoz assemble   \
        --genetic_code 2   \
        --clade Chordata   \
        --outprefix mito_denovo   \
        --thread_number 16   \
        --fq1  ${input_file1}  \
        --fq2 ${input_file2}  \
        --assembler megahit   \
        --kmers_megahit 21 29 39 59 79 99 119 141 \
        --requiring_taxa 9606
    """
}
