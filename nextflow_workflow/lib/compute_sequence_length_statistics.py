import argparse
import sys
import statistics


def get_read_lengths(fastq_file):
    lengths = []
    with open(fastq_file, "r") as f:
        i = 0
        for line in f:
            if i % 4 == 1:  # FASTQ sequence line
                lengths.append(len(line.strip()))
            i += 1

    return lengths


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute sequence length statistics from a FASTQ file."
    )
    parser.add_argument(
        "--fastq_file", type=str, required=True, help="Path to the FASTQ file"
    )
    args = parser.parse_args()
    fq_file = args.fastq_file

    if not fq_file:
        print("Please provide a FASTQ file.")
        sys.exit(1)

    lengths = get_read_lengths(fq_file)
    mean_len = round(statistics.mean(lengths))
    stddev = round(statistics.stdev(lengths))

    # Define both lower and upper length thresholds
    lower_cutoff = mean_len - stddev
    upper_cutoff = mean_len + stddev

    print(f"{lower_cutoff},{upper_cutoff}")
