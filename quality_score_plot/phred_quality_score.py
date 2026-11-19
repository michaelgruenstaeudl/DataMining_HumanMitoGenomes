import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_quality_scores(csv_path: Path | str) -> pd.DataFrame:
    """Load quality scores from one CSV file."""
    return pd.read_csv(csv_path)


def plot_phred_quality_histogram(
    quality_score: pd.DataFrame,
    bins: int = 30,
    output_path: Path | str | None = None,
    show: bool = False,
) -> tuple[plt.Figure, plt.Axes]:
    """Plot Initial_Mean_Phred and Final_Mean_Phred as overlaid histograms.

    Invalid or missing values in either Phred-score column are excluded. Both
    histograms use the same bin edges, making the distributions directly
    comparable.
    """
    if bins < 1:
        raise ValueError("bins must be at least 1")

    columns = ["Initial_Mean_Phred", "Final_Mean_Phred"]
    missing_columns = [column for column in columns if column not in quality_score]
    if missing_columns:
        raise ValueError(f"Missing required columns: {', '.join(missing_columns)}")

    phred_scores = quality_score[columns].apply(pd.to_numeric, errors="coerce").dropna()
    if phred_scores.empty:
        raise ValueError("No valid Phred-score values were found to plot.")

    lower_bound = phred_scores.min().min()
    upper_bound = phred_scores.max().max()
    bin_edges = (
        np.linspace(lower_bound, upper_bound, bins + 1)
        if lower_bound != upper_bound
        else None
    )

    initial_mean = phred_scores["Initial_Mean_Phred"].mean()
    final_mean = phred_scores["Final_Mean_Phred"].mean()
    mean_change = final_mean - initial_mean
    percent_change = mean_change / initial_mean * 100 if initial_mean else float("nan")

    figure, axes = plt.subplots(figsize=(7.2, 4.8))
    axes.hist(
        phred_scores["Initial_Mean_Phred"],
        bins=bin_edges if bin_edges is not None else bins,
        alpha=0.55,
        label=f"Initial (mean = {initial_mean:.2f})",
        color="#0072B2",
        edgecolor="white",
        linewidth=0.5,
    )
    axes.hist(
        phred_scores["Final_Mean_Phred"],
        bins=bin_edges if bin_edges is not None else bins,
        alpha=0.55,
        label=f"Final (mean = {final_mean:.2f})",
        color="#D55E00",
        edgecolor="white",
        linewidth=0.5,
    )
    axes.axvline(
        initial_mean,
        color="#0072B2",
        linestyle="--",
        linewidth=1.5,
        label="Initial overall mean",
    )
    axes.axvline(
        final_mean,
        color="#D55E00",
        linestyle="--",
        linewidth=1.5,
        label="Final overall mean",
    )
    axes.set_title("Change in Mean Phred Quality Score", pad=12)
    axes.set_xlabel("Mean Phred score")
    axes.set_ylabel("Number of samples")
    axes.legend(frameon=False, fontsize=9)
    axes.grid(axis="y", alpha=0.25, linewidth=0.7)
    axes.spines["top"].set_visible(False)
    axes.spines["right"].set_visible(False)
    axes.text(
        0.98,
        0.96,
        f"Overall mean change: {mean_change:+.2f}\n"
        f"({percent_change:+.1f}%; n = {len(phred_scores):,})",
        transform=axes.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "0.8"},
    )
    figure.tight_layout()

    if output_path is not None:
        figure.savefig(output_path, dpi=300, bbox_inches="tight")
    if show:
        plt.show()

    return figure, axes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot initial and final mean Phred scores from one CSV file."
    )
    parser.add_argument("--csv_path", type=Path, help="Input quality-score CSV file")
    parser.add_argument(
        "--output",
        type=Path,
        help="Output image path (defaults to <CSV stem>_phred_comparison.png)",
    )
    parser.add_argument("--bins", type=int, default=30, help="Number of histogram bins")
    parser.add_argument(
        "--show", action="store_true", help="Display the plot interactively"
    )
    arguments = parser.parse_args()

    output_path = arguments.output or arguments.csv_path.with_name(
        f"{arguments.csv_path.stem}_phred_comparison.png"
    )
    plot_phred_quality_histogram(
        load_quality_scores(arguments.csv_path),
        bins=arguments.bins,
        output_path=output_path,
        show=arguments.show,
    )
