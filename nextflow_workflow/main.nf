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

def validate_input_file(String parameter_name, value) {
    if (value == null || value.toString().trim().isEmpty()) {
        error "Parameter '${parameter_name}' must be set to a file."
    }

    def path = file(value)
    if (!path.toFile().isFile()) {
        error "Parameter '${parameter_name}' points to a missing file: ${path}"
    }
    return path
}

def validate_input_directory(String parameter_name, value) {
    if (value == null || value.toString().trim().isEmpty()) {
        error "Parameter '${parameter_name}' must be set to a directory."
    }

    def path = file(value)
    if (!path.toFile().isDirectory()) {
        error "Parameter '${parameter_name}' points to a missing directory: ${path}"
    }
    return path
}

def validate_parameters() {
    if (params.contamination_confidence_threshold < 0 || params.contamination_confidence_threshold > 1) {
        error "Parameter 'contamination_confidence_threshold' must be between 0 and 1."
    }
    if (params.rotation_mismatch_threshold < 0) {
        error "Parameter 'rotation_mismatch_threshold' must be zero or greater."
    }

    def csv_file = validate_input_file('csv_file', params.csv_file)
    def csv_header = csv_file.toFile().withReader { reader -> reader.readLine() }
    def required_columns = ['SRA_Run', 'Nucleotide_AccessionID', 'LibraryLayout']
    def missing_columns = required_columns.findAll { column ->
        !(csv_header ?: '').split(',')*.trim().contains(column)
    }
    if (missing_columns) {
        error "Parameter 'csv_file' is missing required columns: ${missing_columns.join(', ')}"
    }

    [
        reference_fasta: params.reference_fasta,
        config_file: params.config_file,
        rotation_seed_fasta: params.rotation_seed_fasta,
        calculate_sequence_length_threshold_script: params.calculate_sequence_length_threshold_script,
        visualize_file_size_changes_script: params.visualize_file_size_changes_script,
        genome_comparison_plot_script: params.genome_comparison_plot_script,
        trim_mitogenome_duplicate_script: params.trim_mitogenome_duplicate_script
    ].each { parameter_name, value -> validate_input_file(parameter_name, value) }

    [
        fastq_directory: params.fastq_directory,
        gb_directory: params.gb_directory,
        fasta_directory: params.fasta_directory,
        contamination_db: params.contamination_db
    ].each { parameter_name, value -> validate_input_directory(parameter_name, value) }

    if (params.is_calculate_evenness) {
        validate_input_file('merge_pacvr_evenness_figure_script', params.merge_pacvr_evenness_figure_script)
    }

    log.info "Input validation passed for ${csv_file}"
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

    calculate_sequence_length_threshold_script_ch = channel.fromPath(params.calculate_sequence_length_threshold_script)
    rotation_seed_fasta_ch = channel.value(file(params.rotation_seed_fasta))
    trim_mitogenome_duplicate_script_ch = channel.fromPath(params.trim_mitogenome_duplicate_script)
    //Quality control process
    quality_control_workflow(
        sample_reads_ch,
        calculate_sequence_length_threshold_script_ch,
    )
    cutoffs_ch = quality_control_workflow.out.cutoffs_ch

    evenness_calculation_workflow_initial_input_ch = quality_control_workflow.out.repair_read_output_ch
        .join(sample_genbank_ch)
        .join(sample_reference_fasta_ch)

    // quality_control_workflow.out.quality_control_out_ch.view { it -> "QC output: ${it}" }
    size_ch = size_ch.mix(quality_control_workflow.out.repair_read_size_ch)
    size_ch = addSizeTracking(size_ch, quality_control_workflow.out.quality_control_out_ch, "Quality Control Output")

    // Contamination removal and mapping
    contamination_db_channel = channel.fromPath(params.contamination_db)
    contamination_mapping_workflow(
        quality_control_workflow.out.quality_control_out_ch,
        contamination_db_channel,
        contamination_confidence_threshold,
        reference_ch,
    )
    contamination_removal_out_ch = contamination_mapping_workflow.out.contamination_reads_ch
    size_ch = addSizeTracking(size_ch, contamination_removal_out_ch, "Contamination Output")

    mapping_process_output = contamination_mapping_workflow.out.mapping_reads_ch

    size_ch = addSizeTracking(size_ch, mapping_process_output, "Mapping Output")

    if (params.is_calculate_evenness) {
        // Evenness calculation workflow
        //evenness calculation for original reads before cleaning
        evenness_calculation_workflow_initial(
            evenness_calculation_workflow_initial_input_ch,
            'initial_phase',
        )
        initial_evenness_output_ch = evenness_calculation_workflow_initial.out.evenness_output_ch

        // evenness calculation for cleaned reads after mapping
        final_evenness_calculation_input_ch = mapping_process_output
            .join(sample_genbank_ch)
            .join(sample_reference_fasta_ch)
        evenness_calculation_workflow_final(
            final_evenness_calculation_input_ch,
            "mapping_phase",
        )
        cleaned_reads_evenness_output_ch = evenness_calculation_workflow_final.out.evenness_output_ch

    }

    assembly_normalization_workflow(
        mapping_process_output,
        cutoffs_ch,
        config_file_ch,
        trim_mitogenome_duplicate_script_ch,
    )
    //Circular genome rotation process
    rotate_genome_input_ch = assembly_normalization_workflow.out.normalized_fasta_output_ch
        .join(sample_reference_fasta_ch)
        .combine(channel.value(params.rotation_mismatch_threshold))

    rotate_circular_genome_workflow(
        rotate_genome_input_ch,
        rotation_seed_fasta_ch,
    )

    rotated_assembled_ch = rotate_circular_genome_workflow.out.rotated_assembled_ch
    rotated_official_ch = rotate_circular_genome_workflow.out.rotated_official_ch
    rotated_genome_ch = rotate_circular_genome_workflow.out.rotated_genome_ch


    //Haplotype Identification process would go here, after the genome rotation and comparison steps, since it relies on having the final rotated genomes for accurate identification.
    // haplotype_identification_input_ch = rotated_genome_ch
    haplotype_identification_workflow(
        rotated_genome_ch
    )

    report_generation_workflow(
        cutoffs_ch,
        size_ch,
        assembly_normalization_workflow.out.assembled_metadata_info_ch,
        assembly_normalization_workflow.out.normalized_metadata_info_ch,
        rotated_genome_ch,
        rotated_assembled_ch,
        rotated_official_ch,
        haplotype_identification_workflow.out.assembled_haplotype_output_ch,
        haplotype_identification_workflow.out.official_haplotype_output_ch,
        params.is_calculate_evenness,
        params.is_calculate_evenness ? evenness_calculation_workflow_initial.out.evenness_output_ch : channel.empty(),
        params.is_calculate_evenness ? evenness_calculation_workflow_final.out.evenness_output_ch : channel.empty(),
        params.is_calculate_evenness ? channel.fromPath(params.merge_pacvr_evenness_figure_script) : channel.empty(),
        channel.fromPath(params.visualize_file_size_changes_script),
        channel.fromPath(params.genome_comparison_plot_script),
    )

    workflow.onComplete {
        println('✅ Finished all processes!')
    }

    publish:

    // output directories for Quality Control subworkflow
    initial_fastqc_publish_dir = quality_control_workflow.out.initial_fastqc_publish_dir
    fastqc_after_quality_control_publish_dir = quality_control_workflow.out.fastqc_after_quality_control_publish_dir
    repair_output_publish_dir = quality_control_workflow.out.repair_output_publish_dir
    qc_output_publish_dir = quality_control_workflow.out.qc_output_publish_dir

    //Cutoff info CSV
    cutoff_info_csv = report_generation_workflow.out.cutoff_info_csv

    //contamination removal output
    contamination_removal_publish_dir = contamination_mapping_workflow.out.contamination_removal_publish_dir

    //mapping process output
    index_dir = contamination_mapping_workflow.out.index_dir
    mapping_process_publish_dir = contamination_mapping_workflow.out.mapping_process_publish_dir

    //evenness calculation output
    initial_sam_output_publish_dir = params.is_calculate_evenness ? evenness_calculation_workflow_initial.out.sam_output_publish_dir : channel.empty()
    initial_bam_output_publish_dir = params.is_calculate_evenness ? evenness_calculation_workflow_initial.out.bam_output_publish_dir : channel.empty()
    initial_evenness_output_publish_dir = params.is_calculate_evenness ? evenness_calculation_workflow_initial.out.evenness_output_publish_dir : channel.empty()
    final_sam_output_publish_dir = params.is_calculate_evenness ? evenness_calculation_workflow_final.out.sam_output_publish_dir : channel.empty()
    final_bam_output_publish_dir = params.is_calculate_evenness ? evenness_calculation_workflow_final.out.bam_output_publish_dir : channel.empty()
    final_evenness_output_publish_dir = params.is_calculate_evenness ? evenness_calculation_workflow_final.out.evenness_output_publish_dir : channel.empty()
    merged_pacvr_figures_output_publish_dir = report_generation_workflow.out.merged_pacvr_figures_output_publish_dir
    merged_pacvr_figures_files = report_generation_workflow.out.merged_pacvr_figures_files

    //File sizes info CSV
    file_size_csv = report_generation_workflow.out.file_size_csv

    //Genome assembly (mitoZ) output
    mitoz_assembler_publish_dir = assembly_normalization_workflow.out.mitoz_assembler_publish_dir

    //Assembled genome metadata info CSV
    assembled_genome_metadata_csv = report_generation_workflow.out.assembled_genome_metadata_csv

    //normalized genome output
    mito_extend_process_publish_dir = assembly_normalization_workflow.out.mito_extend_process_publish_dir
    normalized_fasta_output_ch = assembly_normalization_workflow.out.normalized_fasta_output_ch

    //normalized genome metadata info CSV
    normalized_genome_metadata_csv = report_generation_workflow.out.normalized_genome_metadata_csv

    //rotated genome output
    rotated_assembled_publish_dir = rotate_circular_genome_workflow.out.rotated_assembled_publish_dir
    rotated_official_publish_dir = rotate_circular_genome_workflow.out.rotated_official_publish_dir

    //rotated genome metadata info CSV
    rotated_genome_metadata_csv = report_generation_workflow.out.rotated_genome_metadata_csv

    //rotated genome CSV for reporting
    rotated_genome_csv_for_reporting = report_generation_workflow.out.rotated_genome_csv_for_reporting

    //haplotype identification output
    assembled_haplotype_output_csv_ch = haplotype_identification_workflow.out.assembled_haplotype_output_csv_ch
    official_haplotype_output_csv_ch = haplotype_identification_workflow.out.official_haplotype_output_csv_ch

    //Haplotype identification metadata info CSV
    assembled_haplotype_info_csv = report_generation_workflow.out.assembled_haplotype_info_csv
    official_haplotype_info_csv = report_generation_workflow.out.official_haplotype_info_csv
}

output {
    // Publish the output directories for Quality Control subworkflow
    initial_fastqc_publish_dir {
        path { meta, _initial_fastqc_dir -> "${meta.id}" }
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

    //cutoff info CSV
    cutoff_info_csv {
        path { _cutoff_info_csv -> "reports" }
    }

    //contamination removal output
    contamination_removal_publish_dir {
        path { meta, _contamination_removal_dir -> "${meta.id}" }
    }

    //mapping process output
    index_dir {
        path { _index_dir -> "." }
    }
    mapping_process_publish_dir {
        path { meta, _mapping_process_publish_dir -> "${meta.id}" }
    }

    //evenness calculation output
    initial_sam_output_publish_dir {
        path { meta, _initial_sam_output_dir -> "${meta.id}" }
    }
    initial_bam_output_publish_dir {
        path { meta, _initial_bam_output_dir -> "${meta.id}" }
    }
    initial_evenness_output_publish_dir {
        path { meta, _initial_evenness_output_dir -> "${meta.id}" }
    }
    final_sam_output_publish_dir {
        path { meta, _final_sam_output_dir -> "${meta.id}" }
    }
    final_bam_output_publish_dir {
        path { meta, _final_bam_output_dir -> "${meta.id}" }
    }
    final_evenness_output_publish_dir {
        path { meta, _final_evenness_output_dir -> "${meta.id}" }
    }
    merged_pacvr_figures_output_publish_dir {
        path { meta, _merged_pacvr_figures_output_dir -> "${meta.id}" }
    }

    //collecting all merged evenness figures in ${outputDir}/combined_coverage_visualization
    merged_pacvr_figures_files {
        path { "." }
    }

    //file sizes info CSV
    file_size_csv {
        path { _file_size_csv -> "reports" }
    }

    //Genome assembly (mitoZ) output
    mitoz_assembler_publish_dir {
        path { meta, _mitoz_assembler_dir -> "${meta.id}" }
    }

    //Assembled genome metadata info CSV
    assembled_genome_metadata_csv {
        path { _assembled_genome_metadata_csv -> "reports" }
    }

    //normalized genome output
    mito_extend_process_publish_dir {
        path { meta, _mito_extend_process_dir -> "${meta.id}" }
    }
    normalized_fasta_output_ch {
        path { meta, _normalized_fasta -> "${meta.id}" }
    }

    //normalized genome metadata info CSV
    normalized_genome_metadata_csv {
        path { _normalized_genome_metadata_csv -> "reports" }
    }

    //rotated genome output
    rotated_assembled_publish_dir {
        path { meta, _rotated_assembled_dir -> "${meta.id}" }
    }
    rotated_official_publish_dir {
        path { meta, _rotated_official_dir -> "${meta.id}" }
    }

    //rotated genome metadata info CSV
    rotated_genome_metadata_csv {
        path { _rotated_genome_metadata_csv -> "reports" }
    }

    //rotated genome CSV for reporting
    rotated_genome_csv_for_reporting {
        path { _rotated_genome_csv_for_reporting -> "reports" }
    }

    //haplotype identification output
    assembled_haplotype_output_csv_ch {
        path { meta, _assembled_haplotype_output_csv_ch -> "${meta.id}" }
    }
    official_haplotype_output_csv_ch {
        path { meta, _official_haplotype_output_csv_ch -> "${meta.id}" }
    }

    //Haplotype identification metadata info CSV
    assembled_haplotype_info_csv {
        path { _assembled_haplotype_info_csv -> "reports" }
    }
    official_haplotype_info_csv {
        path { _official_haplotype_info_csv -> "reports" }
    }
}
