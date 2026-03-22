include { novoplast_process as mito_extend_process } from '../modules/novoplasty_process.nf'
include { trim_mitogenome_duplicate } from '../modules/trim_mitogenome_duplicate.nf'
workflow normalize_complete_genome_length {
    take:
    input_channel
    config_file_ch
    trim_mitogenome_duplicate_script_ch
    parent_output_dir

    main:
    input_channel
        .branch { _sample_id, _sra_read_list, _is_single_end, _read_length, _insert_size, _assembled_fasta, _assembled_topology, assembled_seq_len ->
            too_short: assembled_seq_len as Integer < 15000
            ready: true
        }
        .set { input_ch_branched }

    novoplasty_extend_mito_input = input_ch_branched.too_short.map { sample_id, sra_read_list, is_single_end, read_length, insert_size, assembled_fasta, _assembled_topology, _assembled_seq_len ->
        tuple(sample_id, sra_read_list, is_single_end, read_length, insert_size, assembled_fasta, config_file_ch)
    }

    mito_extend_process(
        novoplasty_extend_mito_input,
        parent_output_dir,
    )
    trim_mitogenome_duplicate_input_ch = mito_extend_process.out.novoplasty_output_ch
        .mix(
            input_ch_branched.ready.map { _sample_id, _sra_read_list, _is_single_end, _read_length, _insert_size, assembled_fasta, _assembled_topology, _assembled_seq_len ->
                tuple(_sample_id, assembled_fasta)
            }
        )
        .combine(trim_mitogenome_duplicate_script_ch)
    // novoplasty_extend_mito_input.view { it -> "Novoplasty extend mito input: ${it}" }

    // mito_extend_process.out.novoplasty_output_ch.view { it -> "Novoplasty output: ${it}" }

    // trim_mitogenome_duplicate_input_ch = mito_extend_process.out.novoplasty_output_ch.combine(trim_mitogenome_duplicate_script_ch)

    // trim_mitogenome_duplicate_input_ch.view { it -> "Trim mitogenome duplicate input: ${it}" }

    trim_mitogenome_duplicate(
        trim_mitogenome_duplicate_input_ch,
        parent_output_dir,
    )

    normalized_fasta_output_ch = trim_mitogenome_duplicate.out.trimmed_fasta_ch

    normalized_fasta_output_ch.view { it -> "Normalized fasta output: ${it}" }

    emit:
    normalized_fasta_output_ch
}
