#!/usr/bin/env nextflow

process repairReads {

    tag "${sample_id}"

    // container 'community.wave.seqera.io/library/bbmap:39.33--60639c9e1473b7a8'
    // container 'oras://community.wave.seqera.io/library/bbmap:39.37--9bf150ff9855d6f6'

    publishDir "results/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(input_file1), path(input_file2), val(encoding_int)
    val phase

    output:
    // val sample_id, emit: sample_id
    tuple val(sample_id), path("repair_read_output_${phase}/${sample_id}_1.fastq", optional: true), path("repair_read_output_${phase}/${sample_id}_2.fastq", optional: true), emit: repaired_reads_output_channel
    path "repair_read_output_${phase}/${sample_id}_singletons.fastq", optional: true

    script:
    def output_dir = "repair_read_output_${phase}"
    """
    mkdir ${output_dir}
    repair.sh \
        in1=${input_file1} \
        in2=${input_file2} \
        out1=${output_dir}/${sample_id}_1.fastq \
        out2=${output_dir}/${sample_id}_2.fastq \
        outs=${output_dir}/${sample_id}_singletons.fastq \
        overwrite=true \
        qin=${encoding_int}
    """
}

workflow {
    input_ch = channel.fromFilePairs("/home/b_thapamagar/BioInformatics/DataMining_HumanMitoGenomes/temp_data2/temp" + '/*_{1,2}.fastq')
        .map { sample_id, read -> tuple(sample_id, read[0], read[1], 33) }
        .view { it -> "Input to repair reads: ${it}" }

    repairReads(input_ch)
    repairReads.out.repaired_reads_output_channel.view { it -> "Output from repair reads: ${it}" }
}
