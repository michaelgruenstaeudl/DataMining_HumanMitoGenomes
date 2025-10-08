// utils.nf or inside main.nf
def trackFileSizes(tuple_input, String processName, String readLabel = "") {
    def sraNumber = tuple_input[0]
    // First element is SRA/sample ID
    def filepath = tuple_input[1]
    def filename = filepath.getName()
    // Second element is File object

    def sizeBytes = filepath.size()
    def sizeKB = (sizeBytes / 1024).round(2)
    def sizeMB = (sizeKB / 1024).round(2)

    return tuple(sraNumber, filename, filepath, processName, readLabel, sizeBytes, sizeKB, sizeMB)
}

def addSizeTracking(size_ch, input_ch, String category) {
    def mapped = input_ch.flatMap { tuple_vals ->
        def sample_id = tuple_vals[0]
        def files = tuple_vals.size() > 1 ? tuple_vals[1..-1] : []

        files
            .withIndex()
            .collect { file, idx ->
                def readLabel = files.size() > 1 ? "Read_${idx + 1}" : "Read_1"
                trackFileSizes(tuple(sample_id, file), category, readLabel)
            }
    }

    return size_ch.concat(mapped)
}
process WriteCSV {

    publishDir "results/reports", mode: 'copy'

    input:
    path csv_tmp

    output:
    path "file_sizes.csv"

    script:
    """
    echo "SRA_Number,File_Name,File_Path,Process_Name,Read_Label,Size_Byte,Size_KB,Size_MB" > file_sizes.csv
    cat ${csv_tmp} >> file_sizes.csv
    """
}
