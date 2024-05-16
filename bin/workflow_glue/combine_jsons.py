#!/usr/bin/env python
"""combine_jsons."""

import json

from ezcharts.components.fastcat import histogram_median
from ezcharts.components.mosdepth import load_mosdepth_summary
import pandas as pd

from .util import wf_parser  # noqa: ABS101
from .report_utils.common import CHROMOSOMES, compute_n50, load_hists  # noqa: ABS101


# Define output metrics names
METRICS = [
    "Contaminated",
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
    "Predicted sex chromosome karyotype",
]


def argparser():
    """Create argument parser."""
    parser = wf_parser("combine_jsons")

    # Required analyses jsons
    parser.add_argument("--jsons", nargs="+", required=False, help="Json files")

    # Required mosdepth inputs
    parser.add_argument(
        "--mosdepth_summary", required=True, help="Mosdepth summary stats"
    )
    parser.add_argument(
        "--mosdepth_thresholds", required=True, help="Mosdepth thresholds"
    )

    # Required bamstats inputs
    parser.add_argument(
        "--bamstats_hists", required=True, help="bamstats reads statistics"
    )
    parser.add_argument(
        "--bamstats_flagstats", required=True, help="bamstats flagstats"
    )

    # Required haplocheck output
    parser.add_argument(
        "--haplocheck",
        required=False,
        help="haplocheck output estimating contamination"
    )

    # additional sample information
    parser.add_argument(
        "--inferred_sex",
        choices=["XX", "XY"],
        help="Genetic sex of sample inferred by the workflow"
    )

    # additional threshold
    parser.add_argument(
        "--contamination_threshold",
        default=0.05,
        help=(
            "Report as contaminated samples with" +
            "haplocheck contamination values above this threshold"
        )
    )

    # Output files
    parser.add_argument("--output", required=True, help="Summary stats file")

    return parser


def hapcheck_contamination(haplocheck):
    """Return contamination level from haplocheck file."""
    hc = pd.read_csv(haplocheck, sep="\t")
    return hc.iat[0, hc.columns.get_loc("Contamination Level")]


def mosdepth_summary_stats(mosdepth_summary, mosdepth_thresholds):
    """Create summary statistics from mosdepth data."""
    _, summary_stats = load_mosdepth_summary(mosdepth_summary)

    # Create new values
    # Average chromosomal coverage
    chrcodes = [f"{key}_region" for key in CHROMOSOMES.keys()]
    chr_mean = summary_stats.loc[summary_stats["chrom"].isin(chrcodes), "mean"].mean()
    # Average total coverage
    total_mean = summary_stats.loc[summary_stats.chrom == "total_region", "mean"].mean()
    # Number of bases over >=N-fold coverage
    # These values starts from col. 5, and change
    # with the required values from the user.
    n_bases_abovex = (
        pd.read_csv(mosdepth_thresholds, sep="\t").iloc[:, 4:].sum().to_dict()
    )
    n_bases_abovex = {k.strip("X"): int(v) for k, v in n_bases_abovex.items()}

    # Create list of metrics and values
    metrics = {
        "Chromosomal depth (mean)": float(chr_mean),
        "Total depth (mean)": float(total_mean),
        "Bases with >=N-fold coverage": n_bases_abovex,
    }
    return metrics


def bamstats_summary_stats(bamstats_hists, bamstats_flagstats):
    """Create summary statistics from bamstats data."""
    # Load length and quality histograms
    len_hist, len_hist_map, len_hist_umap = load_hists(bamstats_hists, "length")
    qual_hist, qual_hist_map, qual_hist_umap = load_hists(bamstats_hists, "quality")
    # Load bamstats flagstats dataset
    flag_stats = pd.read_csv(bamstats_flagstats, sep="\t")

    # Compute the metrics.
    # They are quite a number, so split them here
    # Compute N50 with custom function
    n50 = int(compute_n50(len_hist))
    # Compute yield as sum of read lengths
    read_yield = int(sum(len_hist["start"] * len_hist["count"]))
    # Compute yield of reads with length above size X (X is a multiple of 5K)
    yield_by_base = {}
    for size in range(0, 105000, 5000):
        df = len_hist[len_hist["start"] >= size]
        nx_yield = sum(df["start"] * df["count"])
        yield_by_base[size] = int(nx_yield)
    # Compute number of reads using the read names
    read_number = int(sum(len_hist["count"]))
    # Number of mapped reads
    mapped_reads = int(sum(len_hist_map["count"]))
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
    # Central point of a qual bin
    qual_hist['val'] = 0.5*(qual_hist['start'] + qual_hist['end'])
    # Median read quality
    median_qual = float(histogram_median(qual_hist, x='val'))
    # Median read length
    median_length = int(histogram_median(len_hist))

    # Return values
    metrics = {
        "N50": n50,
        "Yield": read_yield,
        "Yield (reads >=Nbp)": yield_by_base,
        "Read number": read_number,
        "Reads mapped": mapped_reads,
        "Reads mapped (%)": mapped_reads_pct,
        "Primary mappings": primary_aln,
        "Secondary mappings": secondary_aln,
        "Supplementary mappings": suppl_aln,
        "Unmapped reads": unmapped,
        "Median read quality": median_qual,
        "Median read length": median_length,
    }
    return metrics


def main(args):
    """Run entry point."""
    # Create the output dictionary
    out_json = {metric: None for metric in METRICS}

    # Check if the contamination level is >= threshold
    out_json['Contaminated'] = "n/a"
    if args.haplocheck:
        contamination = hapcheck_contamination(args.haplocheck)
        if contamination == 'ND':
            out_json['Contaminated'] = False
        elif contamination == 'NV':
            out_json['Contaminated'] = "n/a"
        else:
            is_contaminated = contamination >= args.contamination_threshold
            out_json['Contaminated'] = True if is_contaminated else False

    # Add bamstats stats
    metrics = bamstats_summary_stats(
        args.bamstats_hists, args.bamstats_flagstats
    )
    out_json.update(metrics)

    # Add mosdepth stats
    metrics = mosdepth_summary_stats(args.mosdepth_summary, args.mosdepth_thresholds)
    out_json.update(metrics)

    # Add results from each analyses json
    if args.jsons:
        for json_file in args.jsons:
            out_json.update(json.load(open(json_file)))

    out_json["Predicted sex chromosome karyotype"] = args.inferred_sex

    # Save output JSON
    json.dump(out_json, open(args.output, "w"))
