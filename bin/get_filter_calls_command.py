#!/usr/bin/env python
"""Construct SV call filtering command."""

import argparse
import sys

threshold_lookup = ['0'] + ['2'] * 10 + ['3'] * 9 + ['5'] * 20 + ['8'] * 100


def import_total_depth(path):
    """Get the average read depth."""
    with open(path, "r") as fh:
        for line in fh:
            if "total" in line:
                return float(line.strip().split("\t")[3])


def parse_arguments():
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--vcf",
        required=True
    )

    parser.add_argument(
        "--target_bedfile",
        help=(
            "Provide the path to a bedfile containing regions"
            " in which to retain SV calls."
        ),
        required=False
    )

    parser.add_argument(
        "--depth_summary",
        help=(
            "Provide the path to a bedfile (e.g. from mosdepth)"
            " containing depth of coverage by region."
        ),
        required=True
    )

    parser.add_argument(
        "--bcftools_threads",
        help=(
            "Number of threads to use for bcftools view (default: 1)."
        ),
        default=1
    )

    parser.add_argument(
        "--min_read_support",
        help=(
            "Set the lower cutoff for read support, "
            "or 'auto' to enable autodetection."
        ),
        required=True,
    )

    parser.add_argument(
        "--min_read_support_limit",
        help=(
            "Set the absolute lower cutoff for read support to"
            " fall back to when using autodetection."
        ),
        required=True,
        type=int
    )

    return parser.parse_args()


def main():
    """Run the entry point."""
    args = parse_arguments()

    # Get min read support filter
    # Todo: Check this
    min_read_support = args.min_read_support_limit
    if args.min_read_support in ['auto']:
        avg_depth = import_total_depth(args.depth_summary)
        avg_depth = min(avg_depth, len(threshold_lookup) - 1) \
            if avg_depth else len(threshold_lookup) - 1
        detected_read_support = int(threshold_lookup[round(avg_depth)])

        if detected_read_support > args.min_read_support_limit:
            min_read_support = detected_read_support

    filter_min_read_support = f'INFO/SUPPORT >= {min_read_support}'

    # Now, make string
    filter_string = f"-i '{filter_min_read_support}'"

    # Add target_bed filter (optional)
    if args.target_bedfile:
        filter_string = f"-T {args.target_bedfile} " + filter_string

    # Print command to stdout
    command = (
        f"bcftools view --threads {args.bcftools_threads} {filter_string} {args.vcf}"
    )
    sys.stdout.write(command)


if __name__ == '__main__':
    main()
