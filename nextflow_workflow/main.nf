#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { calculate_sequence_length_threshold } from './modules/compute_sequence_length_threshold.nf'
include { mapping_process } from './modules/mapping.nf'
include { novoplast_process } from './modules/novoplasty_assembler.nf'
include { contamination_removal } from "./modules/contamination_control.nf"
include { addSizeTracking ; WriteCSV } from './lib/utils.nf'
include { quality_control_workflow } from './sub_workflows/quality_control_workflow.nf'
include { seedExtractionProcess } from './modules/extract_seed.nf'
include {
    evenness_calculation_workflow as evenness_calculation_workflow_initial ;
    evenness_calculation_workflow as evenness_calculation_workflow_final
} from './sub_workflows/evenness_calculation_workflow.nf'
include { merge_pacvr_evenness_figures } from './modules/merge_pacvr_evenness_figure.nf'

workflow {
    all_samples_ch = Channel.fromPath(params.csv_file)
        .splitCsv(header: true, sep: ',')

    initial_channel = all_samples_ch.map { row ->
        tuple(
            row.SRA_Run,
            row.Nucleotide_AccessionID,
            file("${params.fastq_directory}/${row.SRA_Run}_1.fastq"),
            file("${params.fastq_directory}/${row.SRA_Run}_2.fastq"),
            file("${params.gb_directory}/${row.Nucleotide_AccessionID}.gb"),
            file("${params.fasta_directory}/${row.Nucleotide_AccessionID}.fasta"),
        )
    }
    merge_pacvr_evenness_fig_script_ch = channel.fromPath(params.merge_pacvr_evenness_figure_script)

    evenness_calculation_workflow_initial_input_ch = initial_channel.map { sample_id, _accession, fastq1, fastq2, gb_file, fasta_file ->
        tuple(sample_id, fastq1, fastq2, gb_file, fasta_file)
    }
    evenness_calculation_workflow_initial(evenness_calculation_workflow_initial_input_ch, "initial_phase")
    initial_evenness_output_ch = evenness_calculation_workflow_initial.out.evenness_output_ch
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
    // input_ch = Channel.fromFilePairs(params.input_files_directory + '/*_{1,2}.fastq')
    //     .map { sample_id, read -> tuple(sample_id, read[0], read[1]) }
    // .view()

    size_ch = Channel.empty()

    size_ch = addSizeTracking(size_ch, input_ch, "Initial Input Files")

    reference_ch = channel.fromPath(params.reference_fasta)
    seed_mito_ch = channel.fromPath(params.seed_mito)
    config_file_ch = channel.fromPath(params.config_file)
    calculate_sequence_length_threshold_script_ch = channel.fromPath(params.calculate_sequence_length_threshold_script)

    //Quality control process
    quality_control_workflow(input_ch, calculate_sequence_length_threshold_script_ch)
    size_ch = size_ch.mix(quality_control_workflow.out.repair_read_size_ch)
    size_ch = addSizeTracking(size_ch, quality_control_workflow.out.quality_control_out_ch, "Quality Control Output")
    cutoffs_ch = quality_control_workflow.out.cutoffs_ch

    //contamination control process
    contamination_db_channel = channel.fromPath(params.contamination_db)
    contamination_control_input_ch = quality_control_workflow.out.quality_control_out_ch.combine(contamination_db_channel)
    contamination_removal(contamination_control_input_ch)
    size_ch = addSizeTracking(size_ch, contamination_removal.out.contamination_removal_fastq_output, "Contamination Output")


    // Mapping process
    mapping_input_ch = contamination_removal.out.contamination_removal_fastq_output.combine(reference_ch)
    mapping_process(mapping_input_ch)
    size_ch = addSizeTracking(size_ch, mapping_process.out.mapping_process_output, "Mapping Output")


    //Evenness calculation workflow
    final_evenness_calculation_input_ch = mapping_process.out.mapping_process_output.join(gb_channel).join(fasta_channel)
    evenness_calculation_workflow_final(final_evenness_calculation_input_ch, "mapping_phase")
    contamination_removed_evenness_output_ch = evenness_calculation_workflow_final.out.evenness_output_ch

    merge_evenness_figure_input_ch = accession_ch
        .join(initial_evenness_output_ch)
        .join(contamination_removed_evenness_output_ch)
        .combine(merge_pacvr_evenness_fig_script_ch)

    merge_pacvr_evenness_figures(
        merge_evenness_figure_input_ch
    )
    // Write file sizes to CSV
    csv_lines_ch = size_ch.map { tuple ->
        tuple.join(",")
    }
    tmp_csv = csv_lines_ch.collectFile(name: "file_sizes.tmp", newLine: true)

    WriteCSV(tmp_csv)
    size_ch.count().view { "Total tuples in channel: ${it}" }

    // Seed extraction process
    seed_extraction_input_ch = mapping_process.out.mapping_process_output.map { sample_id, mapped_read1, _mapped_read2 ->
        tuple(sample_id, mapped_read1)
    }
    seedExtractionProcess(seed_extraction_input_ch)
    // seedExtractionProcess(quality_control_workflow.out.quality_control_out_ch.map { sample_id, read1, _read2 -> tuple(sample_id, read1) })

    // De novo assembly process
    denovo_assmebly_input_ch = mapping_process.out.mapping_process_output
        .join(
            cutoffs_ch.map { length_cutoffs_ch ->
                def (sample_id, _lower_cutoff, upper_cutoff) = length_cutoffs_ch
                def read_length = upper_cutoff.toInteger() + 1
                def insert_size = upper_cutoff.toInteger() * 2
                tuple(sample_id, read_length, insert_size)
            }
        )
        .combine(seed_mito_ch)
        .combine(config_file_ch)
    // .join(seedExtractionProcess.out.seed_output_channel)


    novoplast_process(denovo_assmebly_input_ch)

    workflow.onComplete {
        println("✅ Finished all processes!")
    }
}
