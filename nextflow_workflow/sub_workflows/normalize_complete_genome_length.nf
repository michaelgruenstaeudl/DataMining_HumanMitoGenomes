include { novoplast_process as mito_extend_process } from '../modules/novoplasty_process.nf'
include { trim_mitogenome_duplicate } from '../modules/trim_mitogenome_duplicate.nf'
workflow normalize_complete_genome_length {
    take:
    input_channel
    parent_output_dir

    main:
    novoplasty_extend_mito_input = input_channel.map { sample_id, reads, is_single_end, read_length, insert_size, assembled_fasta, config_file_ch, _trim_mitogenome_duplicate_script ->
        def reads_list = is_single_end
            ? [reads]
            : reads

        tuple(sample_id, reads_list, is_single_end, read_length, insert_size, assembled_fasta, config_file_ch)
    }
    mito_extend_process(
        novoplasty_extend_mito_input,
        parent_output_dir,
    )

    expand_mitogenome_input_ch = input_channel.map { sample_id, sra_reads_list, is_single_end, assembled_fasta, _trim_mitogenome_duplicate_script_ch ->
        tuple(sample_id, sra_reads_list, is_single_end, assembled_fasta)
    }


    trim_mitogenome_duplicate_input_ch = input_channel.map { sample_id, assembled_fasta, _expand_mitogenome_script_ch, trim_mitogenome_duplicate_script_ch ->
        tuple(sample_id, assembled_fasta, trim_mitogenome_duplicate_script_ch)
    }
    trim_mitogenome_duplicate(
        trim_mitogenome_duplicate_input_ch,
        parent_output_dir,
    )

    normalized_fasta_ch = trim_mitogenome_duplicate.out.trimmed_fasta_ch

    emit:
    normalized_fasta_ch
}
