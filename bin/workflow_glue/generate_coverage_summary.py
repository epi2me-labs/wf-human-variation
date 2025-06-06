#!/usr/bin/env python
"""Coverage summaries."""

import pandas as pd

from .util import wf_parser  # noqa: ABS101


def main(args):
    """Run entry point."""
    df = pd.read_csv(args.mosdepth_threshold, sep='\t')

    # Calculate region sizes
    df["length"] = df["end"] - df["start"]

    # Convert coverage to percentage
    coverage_columns = ["1X", "10X", "15X", "20X", "30X"]
    for col in coverage_columns:
        df[col] = ((df[col] / df["length"]) * 100).round(2)

    # Reorder columns
    df = df[[
        "#chrom", "start", "end", "length", "region", "1X",
        "10X", "15X", "20X", "30X"
    ]]

    averages = pd.read_csv(
        args.mosdepth_average, header=None, sep='\t', names=[
            "#chrom", "start", "end", "region", "avg_coverage"])
    averages['avg_coverage'] = averages['avg_coverage'].round(2)

    merged_df = pd.merge(
        df, averages, on=["#chrom", "start", "end", "region"], how="left")
    merged_df.rename({'#chrom': 'chrom'}, axis=1, inplace=True)

    merged_df.to_csv(args.output, index=False, sep='\t')


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("coverage_summary")
    parser.add_argument(
        "--output", required=True,
        help="Output coverage summary file"
    )
    parser.add_argument(
        "--mosdepth_threshold",
        help="Threshold coverages generated by mosdepth",
    )
    parser.add_argument(
        "--mosdepth_average",
        help="Average coverages generated by mosdepth"
    )
    return parser
