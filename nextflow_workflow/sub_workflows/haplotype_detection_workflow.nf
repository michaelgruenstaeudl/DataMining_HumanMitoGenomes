#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include {
    rotate_genome as rotate_assembled_genome ;
    rotate_genome as rotate_official_genome
} from '../modules/rotate_genome.nf'

include {
    haplotype_identification as assembled_genome_haplotype_detection ;
    haplotype_identification as official_genome_haplotype_detection
} from '../modules/haplotype_identification.nf'

include {
    extract_haplotype_info as extract_assembled_haplotype_info ;
    extract_haplotype_info as extract_official_haplotype_info
} from '../modules/extract_haplotype_info.nf'

workflow haplotype_identification_workflow {
    take:
    haplotype_identification_input_ch
    parent_output_dir

    main:
    asembled_fasta_ch = haplotype_identification_input_ch.map { sample_id, assembled_fasta, _official_fasta ->
        tuple(sample_id, assembled_fasta)
    }

    official_fasta_ch = haplotype_identification_input_ch.map { sample_id, _assembled_fasta, official_fasta ->
        tuple(sample_id, official_fasta)
    }

    assembled_genome_haplotype_detection(asembled_fasta_ch, "assembled_genome", parent_output_dir)
    assembled_haplotype_output_csv_ch = assembled_genome_haplotype_detection.out.haplotype_output_channel
    extract_assembled_haplotype_info(assembled_haplotype_output_csv_ch, parent_output_dir)
    assembled_haplotype_output_ch = extract_assembled_haplotype_info.out.haplotype_info_out_ch

    official_genome_haplotype_detection(official_fasta_ch, "official_genome", parent_output_dir)
    official_haplotype_output_csv_ch = official_genome_haplotype_detection.out.haplotype_output_channel
    extract_official_haplotype_info(official_haplotype_output_csv_ch, parent_output_dir)
    official_haplotype_output_ch = extract_official_haplotype_info.out.haplotype_info_out_ch

    emit:
    assembled_haplotype_output_ch
    official_haplotype_output_ch
}
