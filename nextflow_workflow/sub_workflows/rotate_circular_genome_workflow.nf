#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include {
    rotate_genome as rotate_assembled_genome ;
    rotate_genome as rotate_official_genome
} from '../modules/rotate_genome.nf'
workflow rotate_circular_genome_workflow {
    take:
    rotate_genome_input_ch
    parent_output_dir

    main:
    asembled_fasta_ch = rotate_genome_input_ch.map { sample_id, assembled_fasta, _official_fasta, seed_fasta, n_mismatch ->
        tuple(sample_id, assembled_fasta, seed_fasta, n_mismatch)
    }

    official_fasta_ch = rotate_genome_input_ch.map { sample_id, _assembled_fasta, official_fasta, seed_fasta, n_mismatch ->
        tuple(sample_id, official_fasta, seed_fasta, n_mismatch)
    }

    rotate_assembled_genome(asembled_fasta_ch, parent_output_dir)
    // rotated_assembled_ch = rotate_assembled_genome.out.rotated_output_ch

    // rotated_official_ch = rotate_official_genome.out.rotated_output_ch
    rotate_official_genome(official_fasta_ch, parent_output_dir)
}
