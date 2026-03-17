include { novoplast_process as mito_extend_process } from '../modules/novoplasty_process.nf'
include { trim_mitogenome_duplicate } from '../modules/trim_mitogenome_duplicate.nf'
workflow normalize_complete_genome_length {
    take:
    input_channel
    trim_mitogenome_duplicate_script_ch
    parent_output_dir

    main:
    novoplasty_extend_mito_input = input_channel.map { sample_id, sra_read_list, is_single_end, read_length, insert_size, assembled_fasta, config_file_ch, _trim_mitogenome_duplicate_script ->
        tuple(sample_id, sra_read_list, is_single_end, read_length, insert_size, assembled_fasta, config_file_ch)
    }
    mito_extend_process(
        novoplasty_extend_mito_input,
        parent_output_dir,
    )
    // mito_extend_process.out.novoplasty_output_ch.view { it -> "Novoplasty output: ${it}" }

    trim_mitogenome_duplicate_input_ch = mito_extend_process.out.novoplasty_output_ch.combine(trim_mitogenome_duplicate_script_ch)

    // trim_mitogenome_duplicate_input_ch.view { it -> "Trim mitogenome duplicate input: ${it}" }

    trim_mitogenome_duplicate(
        trim_mitogenome_duplicate_input_ch,
        parent_output_dir,
    )

    normalized_fasta_ch = trim_mitogenome_duplicate.out.trimmed_fasta_ch

    normalized_fasta_ch.view { it -> "Normalized fasta output: ${it}" }

    emit:
    normalized_fasta_ch
}
