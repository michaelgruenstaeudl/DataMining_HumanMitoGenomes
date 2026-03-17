#!/usr/bin/env nextflow

process rotate_genome {

    label "rotate_genome"
    tag "${sample_id}"
    publishDir "${parent_output_dir}/${sample_id}", mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(input_fasta), path(seed_fasta), val(n_mismatch)
    val genome_type
    val parent_output_dir

    output:
    tuple val(sample_id), path(output_file_path), emit: rotated_output_ch
    path "rotated_output/*"

    script:
    output_directory = "rotated_output/${genome_type}"
    output_file_path = "${output_directory}/rotated_${sample_id}.fasta"
    log_file_path = "${output_directory}/rotate.log"
    """
    mkdir -p ${output_directory}
    SEED=\$(grep -v ">" ${seed_fasta})
    rotate -s \$SEED -m ${n_mismatch} ${input_fasta} > ${output_file_path}  2> ${log_file_path}
    
    # Count forward + reverse matches
    MATCH_COUNT=\$(grep -c -E "forward|reverse" ${log_file_path} || echo 0)

    if [ "\$MATCH_COUNT" -gt 1 ]; then
        echo "[WARNING] Multiple matches (\$MATCH_COUNT). Re-rotating using highest position." >> ${log_file_path}

        # Get highest position from all forward/reverse lines
        MAX_POS=\$(grep -E "forward|reverse" ${log_file_path} | awk '{print \$2}' | sort -n | tail -1)
        echo "[INFO] Highest position: \$MAX_POS" >> ${log_file_path}

        # Re-rotate using -x with highest position
        rotate -x \$MAX_POS ${input_fasta}  > ${output_file_path} 2>> ${log_file_path}
        echo "[INFO] Re-rotation complete with -x \$MAX_POS" >> ${log_file_path}
    fi
    """
}
