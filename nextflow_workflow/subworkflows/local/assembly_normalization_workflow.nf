#!/usr/bin/env nextflow

include { mitoz_assembler } from '../../modules/local/mitoz_assembler.nf'
include {
    extract_assembled_genome_metadata_info
    extract_assembled_genome_metadata_info as extract_normalized_genome_metadata_info
} from '../../modules/local/extract_assembled_genome_metadata_info.nf'
include { normalize_complete_genome_length } from './normalize_complete_genome_length.nf'

// Input: mapped tuple(meta, reads), QC cutoffs, and normalization resources.
// Outputs expose assembly, normalized FASTA, metadata, and publish directories.
workflow assembly_normalization_workflow {
    take:
    mapping_reads_ch
    cutoffs_ch
    config_file_ch
    trim_mitogenome_duplicate_script_ch

    main:
    mitoz_assembler(mapping_reads_ch)
    extract_assembled_genome_metadata_info(mitoz_assembler.out.mitoz_assembler_output)

    normalize_fasta_input_ch = mapping_reads_ch
        .join(
            cutoffs_ch.map { meta, _lower_cutoff, upper_cutoff ->
                def read_length = upper_cutoff.toInteger() + 1
                def insert_size = upper_cutoff.toInteger() * 2
                tuple(meta, read_length, insert_size)
            }
        )
        .join(mitoz_assembler.out.mitoz_assembler_output)
        .join(extract_assembled_genome_metadata_info.out.metadata_info_out_ch)

    normalize_complete_genome_length(
        normalize_fasta_input_ch,
        config_file_ch,
        trim_mitogenome_duplicate_script_ch,
    )
    extract_normalized_genome_metadata_info(
        normalize_complete_genome_length.out.normalized_fasta_output_ch
    )

    emit:
    mitoz_assembler_publish_dir = mitoz_assembler.out.mitoz_assembler_publish_dir
    assembled_metadata_info_ch = extract_assembled_genome_metadata_info.out.metadata_info_out_ch
    mito_extend_process_publish_dir = normalize_complete_genome_length.out.mito_extend_process_publish_dir
    normalized_fasta_output_ch = normalize_complete_genome_length.out.normalized_fasta_output_ch
    normalized_metadata_info_ch = extract_normalized_genome_metadata_info.out.metadata_info_out_ch
}
