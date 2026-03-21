import argparse
import os
import sys
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import SeqIO, Align
import pandas as pd

# -------------------------------------------------------------------
# Config
# -------------------------------------------------------------------
RECORDS_PER_FIG = 10  # number of samples to plot per figure (adjust for readability)
PAGE_W = 7  # inches
PAGE_H = 9  # inches


# -------------------------------------------------------------------
# Read FASTA
# -------------------------------------------------------------------
def read_fasta(path):
    record = next(SeqIO.parse(path, "fasta"))
    return str(record.seq).upper(), record.id


# -------------------------------------------------------------------
# Align and classify positions
# -------------------------------------------------------------------
def align_and_classify(assembled_seq, official_seq):
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(assembled_seq, official_seq)[0]
    aligned_assembled = alignment[0]
    aligned_official = alignment[1]
    align_len = len(aligned_assembled)

    match_pos = []
    mismatch_pos = []
    insertion_pos = []  # gap in official   → assembled has EXTRA nucleotide
    deletion_pos = []  # gap in assembled  → assembled is MISSING nucleotide

    for i in range(align_len):
        a = aligned_assembled[i]
        o = aligned_official[i]
        if a == "-" and o != "-":
            deletion_pos.append(
                i
            )  # assembled missing — nucleotide in official but not assembled
        elif o == "-" and a != "-":
            insertion_pos.append(
                i
            )  # assembled extra   — nucleotide in assembled but not official

        elif a == o:
            match_pos.append(i)
        else:
            mismatch_pos.append(i)

    identity_pct = len(match_pos) / align_len * 100

    return {
        "align_len": align_len,
        "match_pos": match_pos,
        "mismatch_pos": mismatch_pos,
        "insertion_pos": insertion_pos,
        "deletion_pos": deletion_pos,
        "identity_pct": identity_pct,
    }


def main(csv_filepath, file_directory):
    # --------------------------------------------------------------------
    # Create output directory if it doesn't exist
    # --------------------------------------------------------------------
    os.makedirs(file_directory, exist_ok=True)
    os.makedirs(f"{file_directory}/reports", exist_ok=True)
    os.makedirs(f"{file_directory}/alignments_plots", exist_ok=True)

    # -------------------------------------------------------------------
    #  Generate figure with one row per sample
    # -------------------------------------------------------------------
    rotated_genome_info_csv = pd.read_csv(csv_filepath)
    n_samples = len(rotated_genome_info_csv)
    n_figures = math.ceil(n_samples / RECORDS_PER_FIG)

    # --------------------------------------------------------------------
    # collect stats for all samples in one dataframe for output
    # --------------------------------------------------------------------
    rotated_genome_info_df = pd.DataFrame(
        columns=[
            "SRA_Number",
            "Assembled_Length",
            "Official_Length",
            "Length_Diff",
            "Identity_Pct",
            "Num_Mismatches",
            "Num_Insertions",
            "Num_Deletions",
        ]
    )

    # --------------------------------------------------------------------
    # Process samples in batches and plot
    # --------------------------------------------------------------------
    for fig_idx in range(n_figures):
        start = fig_idx * RECORDS_PER_FIG
        end = min(start + RECORDS_PER_FIG, n_samples)
        batch = rotated_genome_info_csv[start:end]
        n_rows = len(batch)

        row_h = (
            PAGE_H / RECORDS_PER_FIG
        )  # keeps row density consistent across all pages

        fig, axes = plt.subplots(
            n_rows, 1, figsize=(PAGE_W, n_rows * row_h), squeeze=False
        )

        fig.suptitle(
            "Mitogenome comparison — assembled vs official (all samples)",
            fontsize=15,
        )

        # Process each sample in the batch
        for idx, record in enumerate(batch.itertuples()):
            stat_record = {
                "SRA_Number": record.SRA_Number,
                "Assembled_Length": None,
                "Official_Length": None,
                "Length_Diff": None,
                "Identity_Pct": None,
                "Num_Mismatches": None,
                "Num_Insertions": None,
                "Num_Deletions": None,
            }

            sample_id = record.SRA_Number
            ax = axes[idx][0]

            print(f"[{start + idx + 1}/{n_samples}] Aligning {sample_id} ...")

            assembled_seq, assembled_id = read_fasta(record.rotated_assembled_genome)
            official_seq, official_id = read_fasta(record.rotated_official_genome)

            result = align_and_classify(assembled_seq, official_seq)

            stat_record["Assembled_Length"] = len(assembled_seq)
            stat_record["Official_Length"] = len(official_seq)
            stat_record["Length_Diff"] = abs(len(assembled_seq) - len(official_seq))
            stat_record["Identity_Pct"] = f"{result['identity_pct']:.2f}%"
            stat_record["Num_Mismatches"] = len(result["mismatch_pos"])
            stat_record["Num_Insertions"] = len(result["insertion_pos"])
            stat_record["Num_Deletions"] = len(result["deletion_pos"])
            rotated_genome_info_df.loc[len(rotated_genome_info_df)] = stat_record

            # ------------------------------------------------------------
            # Plot matches, mismatches, insertions, deletions
            # ------------------------------------------------------------
            ax.vlines(
                result["match_pos"], -0.4, 0.4, colors="white", linewidth=0.2, alpha=0.4
            )
            ax.vlines(
                result["mismatch_pos"],
                -0.4,
                0.4,
                colors="firebrick",
                linewidth=0.4,
                alpha=0.9,
            )
            ax.vlines(
                result["insertion_pos"],
                -0.4,
                0.4,
                colors="yellow",
                linewidth=0.4,
                alpha=0.9,
            )
            ax.vlines(
                result["deletion_pos"],
                -0.4,
                0.4,
                colors="blue",
                linewidth=0.4,
                alpha=0.9,
            )

            ax.add_patch(
                plt.Rectangle(
                    (0, -0.4),  # bottom left corner (x, y)
                    result["align_len"],  # width
                    0.8,  # height (-0.4 to 0.4)
                    linewidth=0.8,
                    edgecolor="black",
                    facecolor="none",  # transparent fill
                    clip_on=False,
                    zorder=5,
                )
            )
            ax.set_xlim(0, result["align_len"])
            ax.set_ylim(-1, 1)
            ax.set_yticks([])
            ax.spines[["top", "right", "left"]].set_visible(False)

            # Label on left
            ax.set_ylabel(sample_id, fontsize=10, rotation=0, labelpad=40, va="center")

            # Stats on right
            ax.text(
                1.01,
                0.5,
                f"Original: {len(official_seq):,} bp  \n"
                f"New: {len(assembled_seq):,} bp  \n"
                f"Diff: {abs(len(assembled_seq) - len(official_seq)):,} bp  "
                f"Identity: {result['identity_pct']:.2f}%  \n"
                f"Mismatches: {len(result['mismatch_pos']):,}  \n"
                f"Insertions: {len(result['insertion_pos']):,}  "
                f"Deletions: {len(result['deletion_pos']):,}",
                transform=ax.transAxes,
                fontsize=7,
                va="center",
                color="gray",
            )

            # x label only on last row
            if idx == n_samples - 1:
                ax.set_xlabel("Aligned genome position (bp)", fontsize=9)

        # Legend once at bottom
        legend = [
            mpatches.Patch(facecolor="white", label="Match", edgecolor="black"),
            mpatches.Patch(facecolor="firebrick", label="Mismatch", edgecolor="black"),
            mpatches.Patch(facecolor="yellow", label="Insertion", edgecolor="black"),
            mpatches.Patch(facecolor="blue", label="Deletion", edgecolor="black"),
        ]
        fig.legend(
            handles=legend,
            loc="lower center",
            ncol=4,
            fontsize=9,
            bbox_to_anchor=(0.88, 0.92),
            frameon=True,
            edgecolor="black",
            fancybox=False,
            framealpha=1.0,
            title="Alignment Variants",  # ← legend title
            title_fontsize=10,
        )

        plt.tight_layout()
        plt.savefig(
            f"{file_directory}/alignments_plots/assembly_vs_reference_alignment_part{fig_idx + 1}.png",
            dpi=150,
            bbox_inches="tight",
        )

    # Save stats to CSV
    rotated_genome_info_df.to_csv(
        f"{file_directory}/reports/genome_comparison_stats.csv", index=False
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Genome comparison plotter. It compares assembled mitogenome with official reference and plots the matches, mismatches and gaps along the genome."
    )
    parser.add_argument(
        "--csv_file",
        type=str,
        required=True,
        help="Path to the csv file that contains rotated assembled and rotated official genome paths along with sample IDs",
    )

    parser.add_argument(
        "--directory_name",
        type=str,
        required=True,
        help="directory name to save the plot",
    )

    args = parser.parse_args()
    # fq_file = args.fastq_file
    csv_filepath = args.csv_file
    file_directory = args.directory_name

    if not csv_filepath:
        print("Please provide a CSV file.")
        sys.exit(1)

    main(csv_filepath, file_directory)
