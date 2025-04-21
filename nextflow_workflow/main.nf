#!/usr/bin/env nextflow
include { qualityControl } from './modules/qualityControl.nf'

params.input1 = '/home/b_thapamagar/BioInformatics/NCBIrecordMining/SRA_data/SRR6245232_1.fastq'
params.input2 = '/home/b_thapamagar/BioInformatics/NCBIrecordMining/SRA_data/SRR6245232_2.fastq'

workflow {
    println(" Hello World!  from workflow")
    input_ch1 = Channel.fromPath(params.input1, checkIfExists: true).view()
    input_ch2 = Channel.fromPath(params.input2, checkIfExists: true).view()
    qualityControl(input_ch1, input_ch2)
}
