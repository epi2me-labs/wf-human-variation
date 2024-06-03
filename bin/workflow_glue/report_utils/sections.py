"""Sections for ezcharts report."""
import os

import pandas as pd

import dominate.tags as dom_tags  # noqa: I100,I202

import ezcharts as ezc  # noqa: I202
from ezcharts.components.ezchart import EZChart
from ezcharts.components.fastcat import read_quality_plot
from ezcharts.layout.snippets import Grid, Progress, Stats, Tabs

from .common import THEME  # noqa: ABS101


def get_summary_table(
    hist_df, flagstat_df, n_reads_total, n_bases_total, secondary=False
):
    """Create table with summary statistics.

    :param stats_df: `pd.DataFrame` with bamstats per-read stats
    :param flagstat_df: `pd.DataFrame` with bamstats per-file stats
    :param n_ref_seqs: total number of reference sequences
    :param n_reads_total: total number of reads
    :param n_bases_total: total number of sequenced bases
    :param secondary: _description_, defaults to False
    """
    # get metrics
    if hist_df.empty:
        n_reads = 0
        n_unmapped = 0
        n_bases = 0
        # get percentages
        perc_reads = 0.0
        perc_unmapped = 0.0
        perc_bases = 0.0
    else:
        n_reads = sum(hist_df['count'])
        n_unmapped = flagstat_df["unmapped"].sum()
        n_bases = sum(hist_df['count'] * hist_df['start'])
        # get percentages
        perc_reads = n_reads / n_reads_total * 100
        perc_unmapped = n_unmapped / n_reads * 100
        perc_bases = n_bases / n_bases_total * 100

    # function for progress bar cell showing percentage values
    def percentage_table_cell(value):
        dom_tags.td(
            Progress(
                value_min=0,
                value_max=100,
                value_now=round(value, 1),
                bar_cls="bg-secondary" if secondary else "",
                height=30,
            )
        )

    # create table
    with dom_tags.table(cls="table align-middle"):
        # table header
        with dom_tags.thead():
            dom_tags.th("Metric")
            dom_tags.th("Value")
            dom_tags.th("Percentage")
        # row for reads
        with dom_tags.tr():
            dom_tags.td("Reads")
            dom_tags.td(f"{n_reads:,}")
            percentage_table_cell(perc_reads)
        # unmapped reads
        with dom_tags.tr():
            dom_tags.td("Unmapped reads")
            dom_tags.td(f"{n_unmapped:,}")
            percentage_table_cell(perc_unmapped)
        # bases
        with dom_tags.tr():
            dom_tags.td("Bases")
            dom_tags.td(f"{n_bases:,}")
            percentage_table_cell(perc_bases)


def summary(report, sample_names, hist_df, flagstat_df):
    """Create summary section.

    :param report: report object (`ezcharts.components.reports.labs.LabsReport`)
    :param sample_names: collection of sample names
    :param ref_files: collection of reference file names
    :param ref_seqs: collection of reference sequence names
    :param stats_df: `pd.DataFrame` with bamstats per-read stats
    :param flagstat_df: `pd.DataFrame` with bamstats per-file stats
    """
    with report.add_section("Summary", "Summary"):
        n_reads_total = sum(hist_df['count'])
        n_bases_total = sum(hist_df['count'] * hist_df['start'])
        tabs = Tabs()
        for sample_name in [*sample_names]:
            with tabs.add_tab(
                sample_name if sample_name is not None else "total"
            ):
                # get summary stats for individual sample
                get_summary_table(
                    hist_df,
                    flagstat_df.query(f"sample_name == '{sample_name}'"),
                    n_reads_total,
                    n_bases_total,
                    secondary=True,
                )


def sub_heading(string):
    """Create sub heading in report.

    :param string: string to show in heading
    """
    with dom_tags.div():
        dom_tags.h5(string, {"class": "mb-0 pb-3 pt-3"})


def mapping(report, hists_dir, sample_names):
    """Create alignment stats section.

    :param report: report object (`ezcharts.components.reports.labs.LabsReport`)
    :param stats_df_mapped: `pd.DataFrame` with bamstats per-read stats (aligned reads
        only)
    """
    with report.add_section("Alignment statistics", "Alignment"):
        dom_tags.p(
            """
            This section displays the alignment statistics for the sample processed.
            The left plot shows the mapping accuracy (ranging from 0-100%) vs. the
            number of reads. The right plot shows the depth of sequencing vs. number
            of reads.
            """
        )
        tabs = Tabs()
        # prepare data for cumulative depth plot
        for sample_name in sample_names:
            with tabs.add_tab(sample_name):
                with Grid():
                    with dom_tags.div():
                        # Mapping accuracy
                        acc_hist = pd.read_csv(
                            os.path.join(hists_dir, "accuracy.hist"),
                            sep="\t",
                            names=["start", "end", "count"],
                            dtype={"start": float, "end": float, "count": float},
                        )
                        acc_plt = read_quality_plot(
                            acc_hist, binwidth=1, min_qual=80, max_qual=100)
                        acc_plt.title.text = "Mapping accuracy"
                        acc_plt.xAxis = dict(name="Accuracy (%)", min=80, max=100)
                        EZChart(acc_plt, theme=THEME)
                    with dom_tags.div():
                        # Mapping accuracy
                        cvg_hist = pd.read_csv(
                            os.path.join(hists_dir, "coverage.hist"),
                            sep="\t",
                            names=["start", "end", "count"],
                            dtype={"start": float, "end": float, "count": float},
                        )
                        cvg_plt = read_quality_plot(
                            cvg_hist,
                            binwidth=1,
                            min_qual=0,
                            max_qual=max(cvg_hist['end'])
                        )
                        cvg_plt.title.text = "Coverage"
                        cvg_plt.xAxis = dict(
                            name="Coverage", min=0, max=max(cvg_hist['end']))
                        EZChart(cvg_plt, theme=THEME)


def depths(report, depth_df, sample_name):
    """Create depth section.

    This section contains a plot with depth coverage vs. genome position on the left and
    relative cumulative coverage (i.e. percentage of genome covered to at least a
    certain depth vs. depth) on the right.

    :param report: report object (`ezcharts.components.reports.labs.LabsReport`)
    :param depth_df: `pd.DataFrame` with depth data.
    """
    with report.add_section("Depth of coverage", "Depth"):
        dom_tags.p(
            """
            This section illustrates the depth of coverage of the reference genomes. The
            plot shows coverage vs. genomic position (note that the coordinates on
            the x-axis are the positions along the concatenated reference including all
            reference sequences in the respective reference file).
            """
        )
        # If the dataframe is empty, state it
        if depth_df.empty:
            dom_tags.p(
                """
                The depth file is empty.
                """
            )
        # Otherwise, plot the coverage
        else:
            tabs = Tabs()
            # prepare data for cumulative depth plot
            with tabs.add_tab(sample_name):
                plt = ezc.lineplot(
                    data=depth_df.round(2),
                    x="total_mean_pos",
                    y="depth",
                    hue="chrom",
                )
                plt.title = {"text": "Coverage along reference"}
                plt.xAxis.name = "Position along reference"
                plt.yAxis.name = "Sequencing depth"
                for s in plt.series:
                    s.showSymbol = False
                EZChart(plt, theme=THEME)


def at_a_glance(report, sample_names, values):
    """Generate at-a-glance statistics: read n50 and mean coverage."""
    with report.add_section("At a glance", "Description"):
        # at-a-glace stats + intro
        dom_tags.p(
            """
            This report contains visualisations of statistics that can help in
            understanding the results from the wf-human-variation workflow. Each section
            contains different plots or tables, and in general the results are broken
            down by sample. You can quickly jump to an individual section with the links
            in the header bar.
            """
        )
        # Create tabs for each sample
        tabs = Tabs()
        for sample_name in [*sample_names]:
            with tabs.add_tab(
                sample_name if sample_name is not None else "total"
            ):
                Stats(
                    columns=3,
                    items=[
                        (value, title)
                        for title, value in values.items()
                        ],
                    )
