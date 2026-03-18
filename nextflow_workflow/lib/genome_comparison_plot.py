import argparse
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import SeqIO, Align
import pandas as pd

# -------------------------------------------------------------------
# Config
# -------------------------------------------------------------------
PLOT_OUTPUT = "all_comparisons.png"


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

    match_pos, mismatch_pos, gap_pos = [], [], []

    for i in range(align_len):
        a = aligned_assembled[i]
        o = aligned_official[i]
        if a == "-" or o == "-":
            gap_pos.append(i)
        elif a == o:
            match_pos.append(i)
        else:
            mismatch_pos.append(i)

    identity_pct = len(match_pos) / align_len * 100

    return {
        "align_len": align_len,
        "match_pos": match_pos,
        "mismatch_pos": mismatch_pos,
        "gap_pos": gap_pos,
        "identity_pct": identity_pct,
    }


def main(csv_filepath, file_directory):
    # -------------------------------------------------------------------
    #  Generate figure with one row per sample
    # -------------------------------------------------------------------
    rotated_genome_info_csv = pd.read_csv(csv_filepath)
    n_samples = len(rotated_genome_info_csv)
    rotated_genome_info_df = pd.DataFrame(
        columns=[
            "SRA_Number",
            "Assembled_Length",
            "Official_Length",
            "Length_Diff",
            "Identity_Pct",
            "Num_Mismatches",
            "Num_Gaps",
        ]
    )
    fig, axes = plt.subplots(n_samples, 1, figsize=(18, 1.8 * n_samples), squeeze=False)

    fig.suptitle(
        "Mitogenome comparison — assembled vs official (all samples)",
        fontsize=12,
        y=1.01,
    )

    for record in rotated_genome_info_csv.itertuples():
        stat_record = {
            "SRA_Number": record.SRA_Number,
            "Assembled_Length": None,
            "Official_Length": None,
            "Length_Diff": None,
            "Identity_Pct": None,
            "Num_Mismatches": None,
            "Num_Gaps": None,
        }

        idx = record.Index
        sample_id = record.SRA_Number
        ax = axes[idx][0]

        print(f"[{idx+1}/{n_samples}] Aligning {sample_id} ...")

        assembled_seq, assembled_id = read_fasta(record.rotated_assembled_genome)
        official_seq, official_id = read_fasta(record.rotated_official_genome)

        print(f"  Assembled : {len(assembled_seq):,} bp")
        print(f"  Official  : {len(official_seq):,} bp")
        print(f"  Diff      : {abs(len(assembled_seq) - len(official_seq)):,} bp")

        result = align_and_classify(assembled_seq, official_seq)

        stat_record["Assembled_Length"] = len(assembled_seq)
        stat_record["Official_Length"] = len(official_seq)
        stat_record["Length_Diff"] = abs(len(assembled_seq) - len(official_seq))
        stat_record["Identity_Pct"] = f"{result['identity_pct']:.2f}%"
        stat_record["Num_Mismatches"] = len(result["mismatch_pos"])
        stat_record["Num_Gaps"] = len(result["gap_pos"])
        rotated_genome_info_df.loc[len(rotated_genome_info_df)] = stat_record

        # Plot
        ax.vlines(
            result["match_pos"], -0.4, 0.4, colors="steelblue", linewidth=0.2, alpha=0.4
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
            result["gap_pos"], -0.4, 0.4, colors="orange", linewidth=0.4, alpha=0.9
        )

        ax.set_xlim(0, result["align_len"])
        ax.set_ylim(-1, 1)
        ax.set_yticks([])
        ax.spines[["top", "right", "left"]].set_visible(False)

        # Label on left
        ax.set_ylabel(sample_id, fontsize=8, rotation=0, labelpad=90, va="center")

        # Stats on right
        ax.text(
            1.01,
            0.6,
            f"Assembled: {len(assembled_seq):,} bp  "
            f"Official: {len(official_seq):,} bp  "
            f"Diff: {abs(len(assembled_seq) - len(official_seq)):,} bp",
            transform=ax.transAxes,
            fontsize=7,
            va="center",
            color="gray",
        )
        ax.text(
            1.01,
            0.4,
            f"Identity: {result['identity_pct']:.2f}%  "
            f"Mismatches: {len(result['mismatch_pos']):,}  "
            f"Gaps: {len(result['gap_pos']):,}",
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
        mpatches.Patch(color="steelblue", label="Match"),
        mpatches.Patch(color="firebrick", label="Mismatch"),
        mpatches.Patch(color="orange", label="Gap / indel"),
    ]
    fig.legend(
        handles=legend,
        loc="lower center",
        ncol=3,
        fontsize=9,
        bbox_to_anchor=(0.8, 0.8),
    )

    plt.tight_layout()
    plt.savefig(f"{file_directory}/{PLOT_OUTPUT}", dpi=150, bbox_inches="tight")

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
