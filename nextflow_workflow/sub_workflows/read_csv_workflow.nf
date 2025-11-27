#!/usr/bin/env nextflow

workflow read_csv_workflow {
    take:
    params

    main:
    all_samples_ch = channel.fromPath(params.csv_file)
        .splitCsv(header: true, sep: ',')

    initial_channel = all_samples_ch
        .map { row ->
            def r1 = file("${params.fastq_directory}/${row.SRA_Run}_1.fastq")
            def r2 = file("${params.fastq_directory}/${row.SRA_Run}_2.fastq")
            def gb = file("${params.gb_directory}/${row.Nucleotide_AccessionID}.gb")
            def fa = file("${params.fasta_directory}/${row.Nucleotide_AccessionID}.fasta")

            // Collect missing files
            def missing = []
            if (!r1.exists()) {
                missing << r1
            }
            if (!r2.exists()) {
                missing << r2
            }
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
                file("${params.fastq_directory}/${row.SRA_Run}_1.fastq"),
                file("${params.fastq_directory}/${row.SRA_Run}_2.fastq"),
                file("${params.gb_directory}/${row.Nucleotide_AccessionID}.gb"),
                file("${params.fasta_directory}/${row.Nucleotide_AccessionID}.fasta"),
            )
        }
        .filter { it -> it != null }

    evenness_calculation_workflow_initial_input_ch = initial_channel.map { sample_id, _accession, fastq1, fastq2, gb_file, fasta_file ->
        tuple(sample_id, fastq1, fastq2, gb_file, fasta_file)
    }
    accession_ch = initial_channel.map { sample_id, accession, _fastq1, _fastq2, _gb_file, _fasta_file ->
        tuple(sample_id, accession)
    }
    input_ch = initial_channel.map { sample_id, _accession, fastq1, fastq2, _gb_file, _fasta_file ->
        tuple(sample_id, fastq1, fastq2)
    }
    gb_channel = initial_channel.map { sample_id, _accession, _fastq1, _fastq2, gb_file, _fasta_file ->
        tuple(sample_id, gb_file)
    }
    fasta_channel = initial_channel.map { sample_id, _accession, _fastq1, _fastq2, _gb_file, fasta_file ->
        tuple(sample_id, fasta_file)
    }

    emit:
    evenness_calculation_workflow_initial_input_ch
    accession_ch
    input_ch
    gb_channel
    fasta_channel
}
