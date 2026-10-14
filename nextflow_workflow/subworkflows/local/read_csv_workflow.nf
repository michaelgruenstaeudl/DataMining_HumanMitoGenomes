#!/usr/bin/env nextflow

workflow read_csv_workflow {
    take:
    params

    main:
    all_samples_ch = channel.fromPath(params.csv_file)
        .splitCsv(header: true, sep: ',')



    // Canonical sample input contract:
    // tuple(sample_metadata, reads, genbank, fasta)
    // sample_metadata is a map with String id, String accession, and
    // Boolean single_end. reads is a grouped List<Path>; genbank and fasta
    // are required Path values for the current CSV-driven workflow.
    initial_channel = all_samples_ch
        .map { row ->
            def fastq_list = []
            // Collect missing files
            def missing = []

            def is_single_end = (row.LibraryLayout == "SINGLE")
            if (is_single_end) {
                def r1 = file("${params.fastq_directory}/${row.SRA_Run}*.fastq")
                if (r1.size() == 0) {
                    missing << r1
                }
                fastq_list << r1[0]
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
                [
                    id: row.SRA_Run.toString(),
                    accession: row.Nucleotide_AccessionID.toString(),
                    single_end: is_single_end,
                ],
                fastq_list,
                gb,
                fa,
            )
        }
        .filter { it -> it != null }

    emit:
    sample_input_ch = initial_channel
}
