#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
include {
    rotate_genome as rotate_assembled_genome ;
    rotate_genome as rotate_official_genome
} from '../../modules/local/rotate_genome.nf'

include {
    haplotype_identification as assembled_genome_haplotype_detection ;
    haplotype_identification as official_genome_haplotype_detection
} from '../../modules/local/haplotype_identification.nf'

include {
    extract_haplotype_info as extract_assembled_haplotype_info ;
    extract_haplotype_info as extract_official_haplotype_info
} from '../../modules/local/extract_haplotype_info.nf'

workflow haplotype_identification_workflow {
    take:
    haplotype_identification_input_ch

    main:
    asembled_fasta_ch = haplotype_identification_input_ch.map { meta, assembled_fasta, _official_fasta ->
        tuple(meta, assembled_fasta)
    }

    official_fasta_ch = haplotype_identification_input_ch.map { meta, _assembled_fasta, official_fasta ->
        tuple(meta, official_fasta)
    }

    assembled_genome_haplotype_detection(asembled_fasta_ch, "assembled_genome")
    assembled_haplotype_output_csv_ch = assembled_genome_haplotype_detection.out.haplotype_output_channel
    extract_assembled_haplotype_info(assembled_haplotype_output_csv_ch)
    assembled_haplotype_output_ch = extract_assembled_haplotype_info.out.haplotype_info_out_ch

    official_genome_haplotype_detection(official_fasta_ch, "official_genome")
    official_haplotype_output_csv_ch = official_genome_haplotype_detection.out.haplotype_output_channel
    extract_official_haplotype_info(official_haplotype_output_csv_ch)
    official_haplotype_output_ch = extract_official_haplotype_info.out.haplotype_info_out_ch

    emit:
    assembled_haplotype_output_csv_ch
    official_haplotype_output_csv_ch

    assembled_haplotype_output_ch
    official_haplotype_output_ch
}
