#!/usr/bin/env python
"""combine_jsons."""

import json

from ezcharts.components.mosdepth import load_mosdepth_summary
import pandas as pd

from .util import wf_parser  # noqa: ABS101
from .report_utils.common import CATEGORICAL, CHROMOSOMES, compute_n50  # noqa: ABS101


# Define output metrics names
METRICS = [
    "N50",
    "Yield",
    "Yield (reads >=Nbp)",
    "Read number",
    "Reads mapped",
    "Reads mapped (%)",
    "Primary mappings",
    "Secondary mappings",
    "Supplementary mappings",
    "Unmapped reads",
    "Median read quality",
    "Median read length",
    "Chromosomal depth (mean)",
    "Total depth (mean)",
    "Bases with >=N-fold coverage",
    "SNVs",
    "Indels",
    "Transition/Transversion rate",
    "SV insertions",
    "SV deletions",
    "Other SVs",
]


def argparser():
    """Create argument parser."""
    parser = wf_parser("combine_jsons")

    # Required analyses jsons
    parser.add_argument(
        "--jsons",
        nargs='+',
        required=False,
        help="Json files"
    )

    # Required mosdepth inputs
    parser.add_argument(
        "--mosdepth_summary",
        required=True,
        help="Mosdepth summary stats"
    )
    parser.add_argument(
        "--mosdepth_thresholds",
        required=True,
        help="Mosdepth thresholds"
    )

    # Required bamstats inputs
    parser.add_argument(
        "--bamstats_readstats",
        required=True,
        help="bamstats reads statistics"
    )
    parser.add_argument(
        "--bamstats_flagstats",
        required=True,
        help="bamstats flagstats"
    )

    # Output files
    parser.add_argument("--output", required=True, help="Summary stats file")

    return parser


def mosdepth_summary_stats(mosdepth_summary, mosdepth_thresholds):
    """Create summary statistics from mosdepth data."""
    _, summary_stats = load_mosdepth_summary(mosdepth_summary)

    # Create new values
    # Average chromosomal coverage
    chrcodes = [f"{key}_region" for key in CHROMOSOMES.keys()]
    chr_mean = (
        summary_stats
        .loc[summary_stats['chrom'].isin(chrcodes), "mean"]
        .mean()
    )
    # Average total coverage
    total_mean = (
        summary_stats
        .loc[summary_stats.chrom == "total_region", "mean"]
        .mean()
    )
    # Number of bases over >=N-fold coverage
    # These values starts from col. 5, and change
    # with the required values from the user.
    n_bases_abovex = pd.read_csv(
        mosdepth_thresholds,
        sep="\t"
    ).iloc[:, 4:].sum().to_dict()
    n_bases_abovex = {k.strip('X'): int(v) for k, v in n_bases_abovex.items()}

    # Create list of metrics and values
    metrics = {
        'Chromosomal depth (mean)': float(chr_mean),
        'Total depth (mean)': float(total_mean),
        "Bases with >=N-fold coverage": n_bases_abovex,
    }
    return metrics


def bamstats_summary_stats(bamstats_readstats, bamstats_flagstats):
    """Create summary statistics from bamstats data."""
    read_stats = pd.read_csv(
        bamstats_readstats,
        sep="\t",
        usecols=["name", "ref", "read_length", "mean_quality"],
        dtype={
            "name": str,
            "ref": CATEGORICAL,
            "read_length": int,
            "mean_quality": float
        }
    )
    # Load bamstats flagstats dataset
    flag_stats = pd.read_csv(
        bamstats_flagstats,
        sep="\t"
    )

    # Compute the metrics.
    # They are quite a number, so split them here
    # Compute N50 with custom function
    n50 = int(compute_n50(read_stats.read_length.values))
    # Compute yield as sum of read lengths
    read_yield = int(read_stats.read_length.values.sum())
    # Compute yield of reads with length above size X (X is a multiple of 5K)
    yield_by_base = {}
    for size in range(0, 105000, 5000):
        nx_yield = read_stats[read_stats.read_length >= size].read_length.values.sum()
        yield_by_base[size] = int(nx_yield)
    # Compute number of reads using the read names
    read_number = int(read_stats.name.unique().shape[0])
    # Number of mapped reads
    mapped_reads = int(sum(read_stats.ref != "*"))
    # Number of primary alignments
    primary_aln = int(flag_stats.primary.sum())
    # Number of secondary alignments
    secondary_aln = int(flag_stats.secondary.sum())
    # Number of suppl. alignments
    suppl_aln = int(flag_stats.supplementary.sum())
    # Number of unmapped reads
    unmapped = int(flag_stats.unmapped.sum())
    # Percentage of mapped reads
    mapped_reads_pct = float(round(1 - (unmapped / read_number), 4) * 100)
    # Median read quality
    median_qual = float(read_stats.mean_quality.median())
    # Median read length
    median_length = int(read_stats.read_length.median())

    # Return values
    metrics = {
        'N50': n50,
        'Yield': read_yield,
        'Yield (reads >=Nbp)': yield_by_base,
        'Read number': read_number,
        'Reads mapped': mapped_reads,
        'Reads mapped (%)': mapped_reads_pct,
        'Primary mappings': primary_aln,
        'Secondary mappings': secondary_aln,
        'Supplementary mappings': suppl_aln,
        'Unmapped reads': unmapped,
        'Median read quality': median_qual,
        'Median read length': median_length,
    }
    return metrics


def main(args):
    """Run entry point."""
    # Create the output dictionary
    out_json = {metric: None for metric in METRICS}

    # Add bamstats stats
    metrics = bamstats_summary_stats(
        args.bamstats_readstats,
        args.bamstats_flagstats
    )
    out_json.update(metrics)

    # Add mosdepth stats
    metrics = mosdepth_summary_stats(
        args.mosdepth_summary,
        args.mosdepth_thresholds
    )
    out_json.update(metrics)

    # Add results from each analyses json
    if args.jsons:
        for json_file in args.jsons:
            out_json.update(json.load(open(json_file)))

    # Save output JSON
    json.dump(out_json, open(args.output, 'w'))
