#!/bin/bash
# SRA_NUMBER=DRR058012
#
# csv_file="01_input/common_SRA_in_4_assemblies.csv"

# Access arguments
SRA_NUMBER=$1
ACC_NUMBER=$2
original_figure_file_path=$3
contamination_removed_figure_file_path=$4

# Inputs
original_plot_directory="Data/PACVr_output_Analysis2/original"
normalized_plot_directory="Data/PACVr_output_Analysis2/normalized"
output_directory="combined_coverage_visualization"

latex_file_directory="latex_output"

# temp_file=$(mktemp)
# Skip header and write to temp file
# tail -n +2 "$csv_file" >"$temp_file"
mkdir -p "${latex_file_directory}"

# while IFS=, read -r SRA_NUMBER ACC_NUMBER; do
latex_file_name="${latex_file_directory}/${ACC_NUMBER}.tex"
original_plot_file_path="${original_plot_directory}/${ACC_NUMBER}_PlastomeOnly_original_CoverageViz.pdf"
normalized_plot_file_path="${normalized_plot_directory}/${ACC_NUMBER}_PlastomeOnly_normalized_CoverageViz.pdf"

ESCAPED_ACCESSION=$(echo "$ACC_NUMBER" | sed 's/_/\\_/g')

# Generate LaTeX source
cat >"${latex_file_name}" <<EOF
\documentclass{article}
\usepackage[landscape, top=1in, left=0cm, right=0cm, bottom=0cm,]{geometry} % Landscape + smaller margins
\usepackage{graphicx}
\usepackage{caption}
\usepackage{subcaption} % For side-by-side captions

\pagestyle{empty} % Remove page number

\begin{document}
\begin{center}
    \LARGE \textbf{Accession: ${ESCAPED_ACCESSION}, SRA Number: ${SRA_NUMBER}}
\end{center}
\begin{figure}[!ht]
	\centering
	\begin{subfigure}[t]{0.49\textwidth}
		\centering
		\includegraphics[width=\linewidth, trim=2cm 3cm 4cm 2cm, clip]{${original_figure_file_path}}
		\caption*{\hspace{0.5in}\bfseries Original Reads}
	\end{subfigure}
	\hfill
	\vrule width 1pt\relax
	\hfill
	\begin{subfigure}[t]{0.49\textwidth}
		\centering
		\includegraphics[width=\linewidth, trim=2cm 3cm 4cm 2cm, clip]{${contamination_removed_figure_file_path}}
		\caption*{\hspace{0.5in}\bfseries Contamination removed Reads}
	\end{subfigure}
\end{figure}

\end{document}
EOF

# Compile LaTeX
pdflatex -output-directory=${latex_file_directory} "$latex_file_name"

output_file_name="$output_directory/${ACC_NUMBER}_merged.pdf"
# Rename output
mv "${latex_file_directory}/${ACC_NUMBER}.pdf" "${output_file_name}"

# Clean up auxiliary files
rm -rf "${latex_file_directory:?}"/*

echo "✅ Combined PDF created: $output_file_name"

# done <"$temp_file"

# rm "$temp_file"
