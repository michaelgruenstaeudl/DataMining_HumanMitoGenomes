#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { addSizeTracking } from './modules/local/utils.nf'
include { quality_control_workflow } from './subworkflows/local/quality_control_workflow.nf'
include { read_csv_workflow } from './subworkflows/local/read_csv_workflow.nf'
include {
    evenness_calculation_workflow as evenness_calculation_workflow_initial ;
    evenness_calculation_workflow as evenness_calculation_workflow_final
} from './subworkflows/local/evenness_calculation_workflow.nf'
include { contamination_mapping_workflow } from './subworkflows/local/contamination_mapping_workflow.nf'
include { assembly_normalization_workflow } from './subworkflows/local/assembly_normalization_workflow.nf'
include { rotate_circular_genome_workflow } from './subworkflows/local/rotate_circular_genome_workflow.nf'
include { haplotype_identification_workflow } from './subworkflows/local/haplotype_detection_workflow.nf'
include { report_generation_workflow } from './subworkflows/local/report_generation_workflow.nf'
include { write_csv as write_mean_phred_quality_score_csv } from './modules/local/write_csv.nf'
def validate_input_file(parameter_name: String, value) {
    if (value == null || value.toString().trim().isEmpty()) {
        error("Parameter '${parameter_name}' must be set to a file.")
    }

    def path = file(value)
    if (!path.toFile().isFile()) {
        error("Parameter '${parameter_name}' points to a missing file: ${path}")
    }
    return path
}

def validate_input_directory(parameter_name: String, value) {
    if (value == null || value.toString().trim().isEmpty()) {
        error("Parameter '${parameter_name}' must be set to a directory.")
    }

    def path = file(value)
    if (!path.toFile().isDirectory()) {
        error("Parameter '${parameter_name}' points to a missing directory: ${path}")
    }
    return path
}

def validate_parameters() {
    if (params.contamination_confidence_threshold < 0 || params.contamination_confidence_threshold > 1) {
        error("Parameter 'contamination_confidence_threshold' must be between 0 and 1.")
    }
    if (params.rotation_mismatch_threshold < 0) {
        error("Parameter 'rotation_mismatch_threshold' must be zero or greater.")
    }

    def csv_file = validate_input_file('csv_file', params.csv_file)
    def csv_header = csv_file.toFile().withReader { reader -> reader.readLine() }
    def required_columns = ['SRA_Run', 'Nucleotide_AccessionID', 'LibraryLayout']
    def missing_columns = required_columns.findAll { column ->
        !(csv_header ?: '').split(',')*.trim().contains(column)
    }
    if (missing_columns) {
        error("Parameter 'csv_file' is missing required columns: ${missing_columns.join(', ')}")
    }

    [reference_fasta: params.reference_fasta, config_file: params.config_file, rotation_seed_fasta: params.rotation_seed_fasta, calculate_mean_phred_quality_script: params.calculate_mean_phred_quality_script, calculate_sequence_length_threshold_script: params.calculate_sequence_length_threshold_script, visualize_file_size_changes_script: params.visualize_file_size_changes_script, genome_comparison_plot_script: params.genome_comparison_plot_script, trim_mitogenome_duplicate_script: params.trim_mitogenome_duplicate_script].each { parameter_name, value -> validate_input_file(parameter_name, value) }

    [fastq_directory: params.fastq_directory, gb_directory: params.gb_directory, fasta_directory: params.fasta_directory, contamination_db: params.contamination_db].each { parameter_name, value -> validate_input_directory(parameter_name, value) }

    if (params.is_calculate_evenness) {
        validate_input_file('merge_pacvr_evenness_figure_script', params.merge_pacvr_evenness_figure_script)
    }

    log.info("Input validation passed for ${csv_file}")
}


workflow {

    main:
    if (params.validate_inputs) {
        validate_parameters()
    }

    contamination_confidence_threshold = params.contamination_confidence_threshold
    config_file_ch = file(params.config_file)
    // Call the read_csv_workflow with params
    read_csv_workflow(params)
    sample_input_ch = read_csv_workflow.out.sample_input_ch

    // Compatibility views for process interfaces that still consume
    // sample-specific tuples. The canonical input remains sample_input_ch.
    sample_reads_ch = sample_input_ch.map { sample_metadata, reads, _genbank, _fasta ->
        tuple(sample_metadata, reads)
    }
    sample_genbank_ch = sample_input_ch.map { sample_metadata, _reads, genbank, _fasta ->
        tuple(sample_metadata, genbank)
    }
    sample_reference_fasta_ch = sample_input_ch.map { sample_metadata, _reads, _genbank, fasta ->
        tuple(sample_metadata, fasta)
    }

    size_ch = channel.empty()

    size_ch = addSizeTracking(size_ch, sample_reads_ch, "Initial Input Files")
    // Created once here and passed to mapping workflow, which will then pass it to the mapping processes.
    //This way we avoid creating multiple index files if the same reference is used for multiple samples.
    reference_ch = channel.value(file(params.reference_fasta))

    calculate_mean_phred_quality_script_ch = channel.fromPath(params.calculate_mean_phred_quality_script)
    calculate_sequence_length_threshold_script_ch = channel.fromPath(params.calculate_sequence_length_threshold_script)
    rotation_seed_fasta_ch = channel.value(file(params.rotation_seed_fasta))
    trim_mitogenome_duplicate_script_ch = channel.fromPath(params.trim_mitogenome_duplicate_script)
    //Quality control process
    quality_control_workflow(
        sample_reads_ch,
        calculate_sequence_length_threshold_script_ch,
        calculate_mean_phred_quality_script_ch,
    )
    mean_phred_quality_score_ch = quality_control_workflow.out.mean_phred_quality_score_ch

    mean_phred_quality_score_lines_ch = mean_phred_quality_score_ch.map { meta, initial_phred_score, initial_reads_count, final_phred_score, final_reads_count ->
        [meta.id, initial_phred_score, initial_reads_count, final_phred_score, final_reads_count].join(',')
    }
    mean_phred_quality_score_tmp_csv = mean_phred_quality_score_lines_ch.collectFile(name: 'mean_phred_quality_score.tmp', newLine: true)
    write_mean_phred_quality_score_csv(mean_phred_quality_score_tmp_csv, 'mean_phred_quality_score', 'mean_phred_quality_score.csv')


    workflow.onComplete {
        println('✅ Finished all processes!')
    }

    publish:

    // output directories for Quality Control subworkflow
    initial_fastqc_publish_dir = quality_control_workflow.out.initial_fastqc_publish_dir
    fastqc_after_repair_publish_dir = quality_control_workflow.out.fastqc_process_after_repair_publish_dir
    fastqc_after_quality_control_publish_dir = quality_control_workflow.out.fastqc_after_quality_control_publish_dir
    repair_output_publish_dir = quality_control_workflow.out.repair_output_publish_dir
    qc_output_publish_dir = quality_control_workflow.out.qc_output_publish_dir

    //mean phred quality score CSV
    mean_phred_quality_score_csv = write_mean_phred_quality_score_csv.out.csv_output
}

output {
    // Publish the output directories for Quality Control subworkflow
    initial_fastqc_publish_dir {
        path { meta, _initial_fastqc_dir -> "${meta.id}" }
    }
    fastqc_after_repair_publish_dir {
        path { meta, _fastqc_after_repair_dir -> "${meta.id}" }
    }
    fastqc_after_quality_control_publish_dir {
        path { meta, _fastqc_after_quality_control_dir -> "${meta.id}" }
    }
    repair_output_publish_dir {
        path { meta, _repair_output_dir -> "${meta.id}" }
    }
    qc_output_publish_dir {
        path { meta, _qc_output_publish_dir -> "${meta.id}" }
    }

    //mean phred quality score CSV
    mean_phred_quality_score_csv {
        path { _mean_phred_quality_score_csv -> "reports" }
    }
}
