#!/usr/bin/env nextflow

process calculate_evenness {

    tag "${sample_id}"

    // container 'evolbioinfo/schmutzi:v1.5.6'
    // container "apptainer/pacvr.sif" 
    publishDir "results/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(gb_file), path(bam_file)
    val phase

    output:
    tuple val(sample_id), path("pacvr_${phase}/${sample_id}_${phase}_CoverageViz.pdf"), emit: evenness_output_ch
    path "pacvr_${phase}/*"

    script:
    def output = "pacvr_${phase}/${sample_id}_${phase}_CoverageViz.pdf"

    """
    ## MAPPING ORIGINAL READS
    mkdir -p pacvr_${phase}
    PACVR_LOC=\$(Rscript -e 'cat(find.package("PACVr"))')
    ## VISUALIZATION FOR ORIGINAL READS
    Rscript "\${PACVR_LOC}"/extdata/PACVr_Rscript.R \
        -k ${gb_file} \
        -b ${bam_file} \
        -i 0 \
        -c TRUE \
        -o ${output}

    """
}
