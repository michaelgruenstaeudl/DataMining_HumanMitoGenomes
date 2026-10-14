#!/usr/bin/env nextflow

process calculate_evenness {

    tag "${meta.id}"

    input:
    tuple val(meta), path(gb_file), path(bam_file)
    val phase

    output:
    tuple val(meta), path("pacvr_${phase}/${meta.id}_${phase}_CoverageViz.pdf"), emit: evenness_output_ch
    tuple val(meta), path("pacvr_${phase}"), emit: evenness_output_publish_dir

    script:
    def output = "pacvr_${phase}/${meta.id}_${phase}_CoverageViz.pdf"

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
