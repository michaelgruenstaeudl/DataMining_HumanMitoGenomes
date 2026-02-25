#!/usr/bin/env nextflow

workflow read_csv_workflow {
    take:
    params

    main:
    all_samples_ch = channel.fromPath(params.csv_file)
        .splitCsv(header: true, sep: ',')



    //input channel with tuple of sample_id, accession, fastq1, fastq2, gb_file, fasta_file
    initial_channel = all_samples_ch
        .map { row ->
            def fastq_list = []
            // Collect missing files
            def missing = []

            def is_single_end = (row.LibraryLayout == "SINGLE")
            if (is_single_end) {
                def r1 = file("${params.fastq_directory}/${row.SRA_Run}*.fastq")
                if (!r1.exists()) {
                    missing << r1
                }
                fastq_list << r1
            }
            else {
                def r1 = file("${params.fastq_directory}/${row.SRA_Run}_1.fastq")
                def r2 = file("${params.fastq_directory}/${row.SRA_Run}_2.fastq")
                if (!r1.exists()) {
                    missing << r1
                }
                if (!r2.exists()) {
                    missing << r2
                }
                fastq_list << r1
                fastq_list << r2
            }
            def gb = file("${params.gb_directory}/${row.Nucleotide_AccessionID}.gb")
            def fa = file("${params.fasta_directory}/${row.Nucleotide_AccessionID}.fasta")

            if (!gb.exists()) {
                missing << gb
            }
            if (!fa.exists()) {
                missing << fa
            }

            // If missing, log them and skip row
            if (missing) {
                log.warn("Skipping sample ${row.SRA_Run}. Missing files: ${missing.join(', ')}")
                return null
            }

            return tuple(
                row.SRA_Run,
                row.Nucleotide_AccessionID,
                fastq_list,
                file("${params.gb_directory}/${row.Nucleotide_AccessionID}.gb"),
                file("${params.fasta_directory}/${row.Nucleotide_AccessionID}.fasta"),
                is_single_end,
            )
        }
        .filter { it -> it != null }

    evenness_calculation_workflow_initial_input_ch = initial_channel.map { sample_id, _accession, fastq_list, gb_file, fasta_file, is_single_end ->
        tuple(sample_id, fastq_list, is_single_end, gb_file, fasta_file)
    }
    accession_ch = initial_channel.map { sample_id, accession, _fastq_list, _gb_file, _fasta_file, _is_single_end ->
        tuple(sample_id, accession)
    }
    input_ch = initial_channel.map { sample_id, _accession, fastq_list, _gb_file, _fasta_file, is_single_end ->
        tuple(sample_id, fastq_list, is_single_end)
    }
    gb_channel = initial_channel.map { sample_id, _accession, _fastq_list, gb_file, _fasta_file, _is_single_end ->
        tuple(sample_id, gb_file)
    }
    fasta_channel = initial_channel.map { sample_id, _accession, _fastq_list, _gb_file, fasta_file, _is_single_end ->
        tuple(sample_id, fasta_file)
    }

    emit:
    evenness_calculation_workflow_initial_input_ch
    accession_ch
    input_ch
    gb_channel
    fasta_channel
}
