process write_csv {

    publishDir "${parent_output_dir}/reports", mode: 'copy'

    input:
    path csv_tmp
    val csv_type
    val output_file_name
    val parent_output_dir

    output:
    path (output_file_name), emit: csv_output

    script:
    def csv_header = ""
    if (csv_type == "file_sizes") {
        csv_header = "SRA_Number,File_Name,File_Path,Process_Name,Read_Label,Is_Single_End,Size_Byte,Size_KB,Size_MB"
    }
    else if (csv_type == "assembled_genome_metadata") {
        csv_header = "SRA_Number,Topology,Sequence_Length"
    }
    else if (csv_type == "rotated_genome_info") {
        csv_header = "SRA_Number,rotated_assembled_genome,rotated_official_genome"
    }

    // def file_name = is_file_sizes_csv ? "file_sizes.csv" : "assembled_genome_metadata.csv"
    """
    echo ${csv_header} > ${output_file_name}
    cat ${csv_tmp} >> ${output_file_name}
    """
}
