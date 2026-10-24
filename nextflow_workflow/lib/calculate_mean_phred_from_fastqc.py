#!/usr/bin/env python3

import argparse
import zipfile
from pathlib import Path


def extract_per_sequence_quality(
    zip_path: Path,
) -> tuple[float, int]:
    """Return weighted quality sum and read count from a FastQC ZIP."""

    with zipfile.ZipFile(zip_path) as archive:
        report_names = [
            name for name in archive.namelist() if name.endswith("fastqc_data.txt")
        ]

        if len(report_names) != 1:
            raise ValueError(
                f"Expected one fastqc_data.txt in {zip_path}, "
                f"found {len(report_names)}."
            )

        report = archive.read(report_names[0]).decode("utf-8")

    in_section = False
    weighted_quality_sum = 0.0
    total_reads = 0

    for line in report.splitlines():
        line = line.strip()

        if line.lower().startswith(">>per sequence quality scores"):
            in_section = True
            continue

        if in_section and line == ">>END_MODULE":
            break

        if not in_section or not line or line.startswith("#"):
            continue

        fields = line.split("\t")
        if len(fields) < 2:
            continue

        quality = float(fields[0])
        count = int(float(fields[1]))

        weighted_quality_sum += quality * count
        total_reads += count

    if not in_section:
        raise ValueError(f"Per Sequence Quality Scores section not found in {zip_path}")

    if total_reads == 0:
        raise ValueError(f"No reads found in {zip_path}")

    return weighted_quality_sum, total_reads


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Estimate mean per-read Phred quality from one "
            "single-end or two paired-end FastQC ZIP files."
        )
    )
    parser.add_argument(
        "--fastqc_zip",
        nargs="+",
        type=Path,
        help="One ZIP for single-end or two ZIPs for paired-end data.",
    )
    args = parser.parse_args()

    if len(args.fastqc_zip) not in (1, 2):
        parser.error("Provide one ZIP for single-end or two for paired-end.")

    weighted_sum = 0.0
    total_reads = 0

    for zip_path in args.fastqc_zip:
        if not zip_path.is_file():
            parser.error(f"File not found: {zip_path}")

        quality_sum, read_count = extract_per_sequence_quality(zip_path)
        weighted_sum += quality_sum
        total_reads += read_count

        # print(f"File: {zip_path}")
        # print(f"  Reads: {read_count:,}")

    mean_quality = weighted_sum / total_reads

    # print("\nCombined result")
    # print(f"Total reads: {total_reads:,}")
    print(f"{mean_quality:.4f},{total_reads}")


if __name__ == "__main__":
    main()
