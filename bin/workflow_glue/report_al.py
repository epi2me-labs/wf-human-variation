#!/usr/bin/env python
"""Report using ezcharts."""
from dominate.tags import b, br, p
from ezcharts.components.common import fasta_idx
from ezcharts.components.fastcat import load_bamstats_flagstat
from ezcharts.components.fastcat import SeqSummary
from ezcharts.components.mosdepth import load_mosdepth_regions, load_mosdepth_summary
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Tabs
import pandas as pd

from .report_utils import sections  # noqa: ABS101
from .report_utils import common  # noqa: ABS101

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

    # Load length histograms
    len_hist, len_hist_map, len_hist_umap = common.load_hists(args.hists_dir, 'length')

    # Compute number of reads and read N50
    read_number_str = str(sum(len_hist['count']))
    read_n50 = common.compute_n50(len_hist)

    # get read n50
    n50_str = str(read_n50) + " bp"

    # get total region depth summary
    depth_t, depth_r = load_mosdepth_summary(args.summary_dir)

    # get mean coverage
    logger.info("Compute average cvg...")
    cov = depth_r.loc[depth_r["chrom"] == "total_region", "mean"].values[0]
    cov_str = "{}x".format(cov)

    # read input flagstats data
    flagstat_df = load_bamstats_flagstat(args.flagstat_fn)

    # Define sample name
    sample_names = [args.sample_name]

    # Import depth files when provided, otherwise make an empty df
    # If the dataframe has no columns, it will return a KeyError when
    # subsetting the dataframe. Catch the error and return an empty dataframe.
    try:
        if args.reference_fai:
            depth_df = load_mosdepth_regions(
                args.depths_dir, faidx=faidx, winsize=args.window_size)[[
                    'chrom', 'total_mean_pos', 'depth'
                ]]
        else:
            depth_df = load_mosdepth_regions(args.depths_dir)[[
                'chrom', 'total_mean_pos', 'depth'
            ]]
    except KeyError:
        depth_df = pd.DataFrame()

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
        args.workflow_version
    )

    # If low-cov provided, then display the error
    if args.low_cov:
        with report.add_section("Sample failing", "Fail"):
            p(
                f"""This dataset was not processed by the workflow as it did not
                meet the minimum bam coverage of {args.low_cov}x required.
                """)

    # Add at-a-glance section with intro
    fields = {
        "Sample total reads": read_number_str,
        "Sample read N50": n50_str,
        "Sample mean coverage": cov_str
    }
    sections.at_a_glance(report, sample_names, fields)

    # Add summary table of the input flagstats
    sections.summary(report, sample_names, len_hist, flagstat_df)

    # Combine multiple input files
    with report.add_section("Read statistics", "Stats"):
        p(
            "This section displays the read statistics for the sample processed. A ",
            "description of each plot is as follows:",
            br(), br(), b("Read quality: "),
            "read quality (range cropped to 4-30) vs. the number of reads.", br(),
            b("Read length: "), "read length vs. the number of reads.", br(),
            b("Base yield above read length: "),
            "base yield above a given read length.",
            br(), b("Accuracy: "),
            " mapping accuracy (ranging from 0-100%) vs. the number of reads.",
            br(), b("Coverage: "),
            "coverage (proportion of read spanned by alignment) vs. the number of ",
            "reads."
        )
        tabs = Tabs()
        # prepare data for cumulative depth plot
        for sample_name in sample_names:
            with tabs.add_tab(sample_name):
                SeqSummary(seq_summary=args.hists_dir, alignment_stats=True)

    # Add depth plots
    if not depth_df.empty:
        sections.depths(report, depth_df, args.sample_name)

    # write the report to the output file
    report.write(f"{args.name}")

    logger.info(f"Written report to '{args.name}'.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument(
        "--name",
        help="report name",
    )
    parser.add_argument(
        "--stats_fn",
        help="`bamstats` per-read stats for the sample",
    )
    parser.add_argument(
        "--flagstat_fn",
        help="Directory with `bamstats` per-file stats",
    )
    parser.add_argument(
        "--hists_dir",
        help="Directory with `bamstats` histogram files",
    )
    parser.add_argument(
        "--depths_dir",
        help="Directory with depth files for the sample",
    )
    parser.add_argument(
        "--summary_dir",
        help="Directory with depth summary files for the sample"
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
        type=float,
        help="Define if the QC report should be for low-coverage BAM"
    )
    parser.add_argument(
        "--sample_name",
        type=str,
        help="Name of the sample"
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
    parser.add_argument(
        "--workflow_version", required=True,
        help="Workflow version",
    )
    return parser
