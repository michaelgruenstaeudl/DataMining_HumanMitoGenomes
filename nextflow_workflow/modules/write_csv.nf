process write_csv {

    publishDir "${parent_output_dir}/reports", mode: 'copy'

    input:
    path csv_tmp
    val is_file_sizes_csv
    val output_file_name
    val parent_output_dir

    output:
    path (output_file_name), emit: csv_output

    script:
    def csv_header = is_file_sizes_csv
        ? "SRA_Number,File_Name,File_Path,Process_Name,Read_Label,Size_Byte,Size_KB,Size_MB"
        : "SRA_Number,Topology,Sequence_Length"

    // def file_name = is_file_sizes_csv ? "file_sizes.csv" : "assembled_genome_metadata.csv"
    """
    echo ${csv_header} > ${output_file_name}
    cat ${csv_tmp} >> ${output_file_name}
    """
}
