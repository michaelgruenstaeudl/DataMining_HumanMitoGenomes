import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap


def main(csv_filepath, file_directory):
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

    # Set index
    plot_pivot_df = pivot_df.set_index(["SRA_Number", "Read_Label"])

    # Prepare data
    categories = plot_pivot_df.columns
    index_labels = plot_pivot_df.index.tolist()  # .sort()
    values = plot_pivot_df.values

    # Positions for bars
    y_pos = np.arange(len(index_labels))

    # Plot
    fig, ax = plt.subplots(figsize=(20, 13))
    # custom_colors = ["red", "black", "green", "yellow", "blue"]
    # custom_colors= ["#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3"]

    custom_colors = ["#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854"]

    my_cmap = ListedColormap(custom_colors)
    colors = my_cmap(np.linspace(0, 1, len(categories)))

    # Draw stacked bars manually
    left = np.zeros(len(index_labels))
    for i, cat in enumerate(categories):
        ax.barh(y_pos, values[:, i], height=0.3, left=left, color=colors[i], label=cat)
        left += values[:, i]

    # Customize
    # ax.invert_yaxis()
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"{sra} | {read}" for sra, read in index_labels])
    # ax.set_xlabel('Value')

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
    # plt.show()


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
