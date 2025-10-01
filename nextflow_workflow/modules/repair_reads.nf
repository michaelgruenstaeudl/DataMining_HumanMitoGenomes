#!/usr/bin/env nextflow

process repairReads {

    tag "${sample_id}"

    publishDir "results/${sample_id}", mode: 'copy', overwrite: true
    container 'community.wave.seqera.io/library/bbmap:39.33--60639c9e1473b7a8'

    input:
    tuple val(sample_id), path(input_file1), path(input_file2), val(encoding_int)

    output:
    // val sample_id, emit: sample_id
    tuple val(sample_id), path("repair_read_output/${sample_id}_1.fastq", optional: true), path("repair_read_output/${sample_id}_2.fastq", optional: true), emit: repaired_reads_output_channel
    path "repair_read_output/${sample_id}_singletons.fastq", optional: true

    script:
    """
    mkdir repair_read_output
    repair.sh in1=${input_file1} in2=${input_file2} out1=repair_read_output/${sample_id}_1.fastq out2=repair_read_output/${sample_id}_2.fastq outs=repair_read_output/${sample_id}_singletons.fastq overwrite=true qin=${encoding_int}
    """
}

workflow {
    input_ch = Channel.fromFilePairs("/home/b_thapamagar/BioInformatics/DataMining_HumanMitoGenomes/temp_data2/temp" + '/*_{1,2}.fastq')
        .map { sample_id, read -> tuple(sample_id, read[0], read[1], 33) }
        .view { "Input to repair reads: ${it}" }

    repairReads(input_ch)
    repairReads.out.repaired_reads_output_channel.view { "Output from repair reads: ${it}" }
}
