import argparse
import random
import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D


def calculate_file_size_reduced_after_each_process(file_size_df):
    file_size_df = pd.read_csv(csv_filepath)
    process_order = [
        "Initial Input Files",
        "Repair Read Output",
        "Quality Control Output",
        "Contamination Output",
        "Mapping Output",
    ]
    # Create a ranking dictionary
    order_map = {name: i for i, name in enumerate(process_order)}

    # Add a temporary sort key column
    file_size_df["_process_order"] = file_size_df["Process_Name"].map(order_map)

    ordered_file_size_df = file_size_df.sort_values(
        by=["SRA_Number", "Read_Label", "_process_order"], ascending=[True, True, True]
    ).drop(columns=["_process_order"])

    ordered_file_size_df["diff"] = (
        ordered_file_size_df.groupby(["SRA_Number", "Read_Label"])["Size_MB"]
        .diff(-1)
        .fillna(ordered_file_size_df["Size_MB"])
    )
    ordered_file_size_df["percentage_diff"] = (
        (
            ordered_file_size_df["diff"]
            / ordered_file_size_df.groupby(["SRA_Number", "Read_Label"])[
                "diff"
            ].transform("sum")
        )
        * 100
    ).round(2)

    ordered_file_size_df["diff_after_process"] = (
        ordered_file_size_df.groupby(["SRA_Number", "Read_Label"])["Process_Name"]
        .shift(-1)  # get next row’s process_name within each group
        .fillna("Final remaining size for reassembly")  # if last row
    )
    ordered_file_size_df

    pivot_df = ordered_file_size_df.pivot(
        index=["SRA_Number", "Read_Label"],
        columns="diff_after_process",
        values="percentage_diff",
    )
    pivot_df = pivot_df.rename(
        columns={
            "Repair Read Output": "Size reduced after repair process",
            "Quality Control Output": "Size reduced after Quality Control poccess",
            "Contamination Output": "Size reduced after Contamination Removal poccess",
            "Mapping Output": "Size reduced after Mapping poccess",
            "Remaining file Size for reassembly": "Final remaining size for reassembly",
        }
    )
    pivoted_difference_order = [
        "Size reduced after repair process",
        "Size reduced after Quality Control poccess",
        "Size reduced after Contamination Removal poccess",
        "Size reduced after Mapping poccess",
        "Final remaining size for reassembly",
    ]

    existing_columns = [
        col for col in pivoted_difference_order if col in pivot_df.columns
    ]
    pivot_df = pivot_df[existing_columns]

    pivot_df = pivot_df.reset_index()
    return pivot_df


def plot_stacked_bar_chart(pivot_df, file_directory):
    plot_pivot_df = pivot_df.set_index(["SRA_Number", "Read_Label"]).fillna(0)

    # Prepare data
    categories = plot_pivot_df.columns
    index_labels = plot_pivot_df.index.tolist()  # .sort()
    values = plot_pivot_df.values

    # Positions for bars
    y_pos = np.arange(len(index_labels))

    num_rows = len(index_labels) / 2  # Adjust divisor to control height
    height_per_row = 0.5  # Adjust this value as needed

    fig_height = num_rows * height_per_row

    # Plot
    fig, ax = plt.subplots(figsize=(20, fig_height))
    # custom_colors = ["red", "black", "green", "yellow", "blue"]
    # custom_colors= ["#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3"]

    custom_colors = ["#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854"]

    my_cmap = ListedColormap(custom_colors)
    colors = my_cmap(np.linspace(0, 1, len(categories)))

    # Draw stacked bars manually
    left = np.zeros(len(index_labels))
    for i, cat in enumerate(categories):
        ax.barh(
            y_pos,
            values[:, i],
            height=height_per_row,
            left=left,
            color=colors[i],
            label=cat,
            edgecolor="black",
        )
        left += values[:, i]

    # Customize
    # ax.invert_yaxis()
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"{sra} | {read}" for sra, read in index_labels])
    # ax.set_xlabel('Value')
    ax.set_ymargin(0.001)
    ax.invert_yaxis()
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()
    plt.xlabel("Percentage Difference")
    plt.title("Stacked Bar Plot for each (SRA_Number, Read_Label)")
    # plt.xticks(rotation=45, ha="right")
    plt.legend(title="Size in Percentage", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(
        f"{file_directory}/file_size_percentage_difference_stacked_bar_plot.png",
        dpi=300,
    )


def plot_aggregate_stacked_bar_chart(pivot_df, file_directory):
    aggregated_df = pivot_df.groupby("SRA_Number").mean(numeric_only=True).fillna(0)

    # Prepare data
    categories = aggregated_df.columns
    index_labels = aggregated_df.index.tolist()  # .sort()
    values = aggregated_df.values

    # Positions for bars
    y_pos = np.arange(len(index_labels))

    num_rows = len(index_labels) / 2  # Adjust divisor to control height
    height_per_row = 0.5  # Adjust this value as needed

    fig_height = num_rows * height_per_row

    # Plot
    fig, ax = plt.subplots(figsize=(20, fig_height))

    custom_colors = ["#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854"]

    my_cmap = ListedColormap(custom_colors)
    colors = my_cmap(np.linspace(0, 1, len(categories)))

    # Draw stacked bars manually
    left = np.zeros(len(index_labels))
    for i, cat in enumerate(categories):
        ax.barh(
            y_pos,
            values[:, i],
            height=height_per_row,
            left=left,
            color=colors[i],
            label=cat,
            edgecolor="black",
        )
        left += values[:, i]

    # Customize
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"{sra}" for sra in index_labels])
    ax.set_ymargin(0.001)
    ax.invert_yaxis()
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()
    plt.xlabel("Percentage Difference")
    plt.title("Stacked Bar Plot for each SRA_Number")
    plt.legend(title="Size in Percentage", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(
        f"{file_directory}/aggregate_file_size_percentage_difference_stacked_bar_plot.png",
        dpi=300,
    )


def plot_file_size_changes(csv_filepath, file_directory):
    file_size_df = pd.read_csv(csv_filepath)
    aggdf = (
        file_size_df.groupby(["SRA_Number", "Process_Name"])["Size_MB"]
        .mean()
        .reset_index()
    )

    pivot_aggdf = (
        aggdf.pivot(
            index=["SRA_Number"],
            columns="Process_Name",
            values="Size_MB",
        )
        .sort_index()
        .round(3)
    )

    n_records = len(pivot_aggdf)
    subplots_per_row = 15
    if n_records <= 50:
        subplots_per_row = 5
    elif n_records <= 100:
        subplots_per_row = 10
    else:
        subplots_per_row = 15

    n_rows = math.ceil(n_records / subplots_per_row)

    custom_colors = ["#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854"]

    desired__process_order = [
        "Initial Input Files",
        "Repair Read Output",
        "Quality Control Output",
        "Contamination Output",
        "Mapping Output",
    ]

    color_dict = dict(zip(desired__process_order, custom_colors))

    value_cols = [c for c in desired__process_order if c in pivot_aggdf.columns]
    fig = plt.figure(figsize=(150, 200))

    # Create GridSpec with controlled spacing
    gs = gridspec.GridSpec(
        n_rows,
        subplots_per_row,
        figure=fig,
        wspace=0.3,  # horizontal gap between subplots
        hspace=0.5,  # vertical gap between subplots
        left=0.01,
        right=0.99,
        bottom=0.01,
        top=0.99,
    )

    axs = []
    for i in range(n_rows):
        for j in range(subplots_per_row):
            ax = fig.add_subplot(gs[i, j])
            axs.append(ax)

    for ax, (_, row) in zip(axs, pivot_aggdf.iterrows()):
        values = row[value_cols]
        values = values.dropna()

        # Match colors to remaining columns
        bar_colors = [custom_colors[value_cols.index(col)] for col in values.index]

        ax.bar(
            values.index,
            values.values,
            color=bar_colors,
            width=0.3,
            edgecolor="black",
            linewidth=0.5,
        )

        # --- LINE PLOT ---
        ax.plot(
            values.index,  # same x-axis
            values.values,  # y values
            color="black",  # line color
            marker="o",  # circle markers
            linewidth=2,  # line thickness
            markersize=6,  # marker size
        )

        ax.set_xticks([])
        ax.set_yticks([])

        ax.set_title(str(row.name), fontsize=75)
        # ax.set_ylabel([])
        ax.grid(axis="y", linestyle="--", alpha=0.5)

    for ax in axs[n_records:]:
        ax.axis("off")

    # -------------------------
    # Highlight the example subplot
    # -------------------------
    example_idx = random.randint(0, len(axs) - 1)
    highlight_ax = axs[example_idx]
    rect = Rectangle(
        (0, 0),
        1,
        1,
        linewidth=20,
        edgecolor="red",
        facecolor="none",
        linestyle="--",
        transform=highlight_ax.transAxes,
        zorder=10,
    )
    highlight_ax.add_patch(rect)

    #### Add one large subplot for the example record
    example_ax = fig.add_axes([1.05, 0.4, 0.5, 0.3])  # large subplot at bottom

    example_record = pivot_aggdf.iloc[example_idx].reindex(desired__process_order)

    example_ax.bar(
        example_record.index,
        example_record.values,
        color=bar_colors,
        width=0.3,
        edgecolor="black",
        linewidth=2,
    )

    example_ax.plot(
        example_record.index,  # same x-axis
        example_record.values,  # y values
        color="black",  # line color
        marker="o",  # circle markers
        linewidth=2,  # line thickness
        markersize=6,  # marker size
    )

    # Set ticks safely
    example_ax.set_xticks(range(len(example_record)))
    example_ax.set_xticklabels(example_record.index, rotation=30, fontsize=100)

    example_ax.set_yticks(example_ax.get_yticks())
    example_ax.set_yticklabels(
        [f"{y:.1f} MB" for y in example_ax.get_yticks()], fontsize=100
    )

    example_ax.set_title(f"Example Record: {example_record.name}", fontsize=150)
    example_ax.set_ylabel("File Size in MB", fontsize=100)
    example_ax.grid(axis="y", linestyle="--", alpha=0.5)

    # -------------------------
    # Optional: connecting line
    # -------------------------
    fig.canvas.draw()  # update coordinates

    bbox_small = highlight_ax.get_window_extent().transformed(
        fig.transFigure.inverted()
    )

    bbox_example = example_ax.get_window_extent().transformed(
        fig.transFigure.inverted()
    )

    top_line = Line2D(
        [bbox_small.x1, bbox_example.x0],  # right edge → left edge
        [bbox_small.y1, bbox_example.y1],  # top → top
        transform=fig.transFigure,
        color="red",
        linewidth=5,
        linestyle="--",
        zorder=20,
    )

    bottom_line = Line2D(
        [bbox_small.x1, bbox_example.x0],  # right edge → left edge
        [bbox_small.y0, bbox_example.y0],  # bottom → bottom
        transform=fig.transFigure,
        color="red",
        linewidth=5,
        linestyle="--",
        zorder=20,
    )

    fig.add_artist(bottom_line)

    fig.add_artist(top_line)

    # -------------------------
    # Add legend
    # -------------------------
    # Create handles manually

    legend_handles = [
        Patch(color=color_dict[name], label=name) for name in desired__process_order
    ]

    # Place legend at top right of figure (outside main grid)
    fig.legend(
        handles=legend_handles,
        title="Process Name",
        loc="upper right",
        bbox_to_anchor=(1.3, 0.99),
        fontsize=150,
        title_fontsize=200,
        frameon=True,
    )

    # -------------------------
    # Figure caption at bottom
    # -------------------------
    fig.text(
        0.5,
        -0.01,
        "Fig: File Size across Nextflow Processing Steps",
        ha="center",
        fontsize=200,
        fontstyle="italic",
    )

    plt.savefig(
        f"{file_directory}/aggregate_file_size_reduction_viz.png",
        bbox_inches="tight",
        dpi=50,
    )


def main(csv_filepath, file_directory):
    pivot_df = calculate_file_size_reduced_after_each_process(csv_filepath)
    # Set index
    plot_stacked_bar_chart(pivot_df, file_directory)

    plot_aggregate_stacked_bar_chart(pivot_df, file_directory)
    plot_file_size_changes(csv_filepath, file_directory)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Visualize file sizes changes over each processes."
    )
    parser.add_argument(
        "--csv_file",
        type=str,
        required=True,
        help="Path to the csv file that contains file sizes information",
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
