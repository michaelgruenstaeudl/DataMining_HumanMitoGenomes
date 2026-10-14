#!/usr/bin/env nextflow

def trackFileSizes(meta, filepath, processName, readLabel) {
    // println("Tracking file size for sample ${meta.id} in process ${processName} for ${readLabel} with file: ${filepath}")
    def sraNumber = meta.id
    def filename = filepath.getName()
    // Second element is File object

    def sizeBytes = filepath.size()
    def sizeKB = (sizeBytes / 1024).round(2)
    def sizeMB = (sizeKB / 1024).round(2)

    return tuple(sraNumber, filename, filepath, processName, readLabel, meta.single_end, sizeBytes, sizeKB, sizeMB)
}

def addSizeTracking(size_ch, input_ch, category) {
    // input_ch.view { it -> "Input to addSizeTracking for ${category}: ${it}" }
    def mapped = input_ch.flatMap { meta, files ->
        // println("Processing sample ${meta.id} in category ${category} with files: ${files}")
        files
            .withIndex()
            .collect { file, idx ->
                // println("Tracking file size for sample ${meta.id} in category ${category} for file: ${file}")
                def readLabel = files.size() > 1 ? "Read_${idx + 1}" : "Read_1"
                trackFileSizes(meta, file, category, readLabel)
            }
    }

    return size_ch.concat(mapped)
}
