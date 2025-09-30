// utils.nf or inside main.nf
def trackFileSizes(tuple_inpput, String processName, String readLabel = "") {
    def sraNumber = tuple_inpput[0]
    // First element is SRA/sample ID
    def file = tuple_inpput[1]
    // Second element is File object

    def sizeBytes = file.size()
    def sizeKB = (sizeBytes / 1024).round(2)
    def sizeMB = (sizeKB / 1024).round(2)

    return tuple(sraNumber, file, processName, readLabel, sizeBytes, sizeKB, sizeMB)
}

process WriteCSV {

    publishDir "results/reports", mode: 'copy'

    input:
    path csv_tmp

    output:
    path "file_sizes.csv"

    script:
    """
    echo "SRA_Number,File_Name,Process_Name,Read_Label,Size_Byte,Size_KB,Size_MB" > file_sizes.csv
    cat ${csv_tmp} >> file_sizes.csv
    """
}
