#!/usr/bin/env python
"""Report using ezcharts."""
from dominate.tags import p
from ezcharts.components.common import fasta_idx
from ezcharts.components.fastcat import load_bamstats_flagstat, load_stats
from ezcharts.components.fastcat import SeqSummary
from ezcharts.components.mosdepth import load_mosdepth_regions
from ezcharts.components.reports import labs
import pandas as pd

from .report_utils import sections  # noqa: ABS101

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run entry point."""
    logger = get_named_logger("report")

    # Import fai file if specified.
    # CW-2585: if no fai is specified, the depth coverage magnify the region(s)
    # coverage only.
    if args.reference_fai:
        faidx = fasta_idx(args.reference_fai)
        if faidx.empty:
            raise pd.errors.EmptyDataError(f'{args.reference_fai}')

    # read input stats data
    stats_df = load_stats(args.stats_dir, format='bamstats')

    # read input flagstats data
    flagstat_df = load_bamstats_flagstat(args.flagstat_dir)

    # Define categories
    sample_names = stats_df["sample_name"].cat.categories
    if sample_names != sample_names:
        raise ValueError('Sample names in the two stats file do not match')

    # Import depth files when provided, otherwise make an empty df
    if args.reference_fai:
        depth_df = load_mosdepth_regions(
            args.depths_dir, faidx=faidx, winsize=args.window_size)
    else:
        depth_df = load_mosdepth_regions(args.depths_dir)

    # create the report
    if args.low_cov:
        report_name = f"{args.name} reads QC report - failing"
    else:
        report_name = f"{args.name} reads QC report"
    report = labs.LabsReport(
        report_name,
        args.name,
        args.params,
        args.versions,
    )

    # If low-cov provided, then display the error
    if args.low_cov:
        with report.add_section("Sample failing", "Fail"):
            p(
                f"""This dataset was not processed by the workflow as it did not
                meet the minimum bam coverage of {args.low_cov}x required.
                """)

    # Add summary table of the input flagstats
    sections.summary(report, sample_names, stats_df, flagstat_df)

    # Combine multiple input files
    with report.add_section("Read statistics", "Stats"):
        SeqSummary(f"{args.stats_dir}/")

    # extract the mapped reads and some other metrics used in the report sections
    stats_df_mapped = stats_df.query('ref != "*"')
    sections.mapping(report, stats_df_mapped)

    # Add depth plots
    if not depth_df.empty:
        depth_df['sample_name'] = depth_df.filename.str.split('.').str[0]
        sections.depths(report, depth_df)

    # write the report to the output file
    report.write(f"{args.name}-alignment-report.html")

    logger.info(f"Written report to '{args.name}-alignment-report.html'.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument(
        "--name",
        help="report name",
    )
    parser.add_argument(
        "--stats_dir",
        help="directory with `bamstats` per-read stats for the sample",
    )
    parser.add_argument(
        "--flagstat_dir",
        help="directory with `bamstats` per-file stats",
    )
    parser.add_argument(
        "--depths_dir",
        help="directory with depth files for the sample",
    )
    parser.add_argument(
        "--reference_fai",
        help="Reference fai index with sequences lengths",
    )
    parser.add_argument(
        "--window_size",
        default=25000,
        type=int,
        help="Size of windows for the depth plot",
    )
    parser.add_argument(
        "--low_cov",
        type=int,
        help="define if the QC report should be for low-cov bam"
    )
    parser.add_argument(
        "--params",
        default=None,
        help="CSV file with workflow parameters",
    )
    parser.add_argument(
        "--versions",
        help="CSV file with software versions",
    )
    return parser
