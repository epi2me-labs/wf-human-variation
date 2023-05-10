"""Sections for ezcharts report."""
import pandas as pd

import dominate.tags as dom_tags  # noqa: I100,I202

import ezcharts as ezc  # noqa: I202
from ezcharts.components.ezchart import EZChart
from ezcharts.layout.snippets import Grid, Progress, Tabs

from .common import CHROMOSOMES, THEME  # noqa: ABS101


def histogram_with_mean_and_median(
    series,
    title=None,
    x_axis_name=None,
    y_axis_name=None,
    bins=100,
    round_digits=1,
):
    """Create ezcharts histogram showing the mean and median underneath the plot title.

    :param series: `pd.Series` with data to plot
    :param title: plot title, defaults to None
    :param x_axis_name: x axis label, defaults to None
    :param y_axis_name: y axis label, defaults to None
    :param bins: number of bins, defaults to 100
    :param round_digits: number of decimals to round the mean and median values to,
        defaults to 1
    :raises ValueError: Raise error if `series` is not a `pd.Series`
    :return: the histogram (`ezcharts.plots.Plot`)
    """
    if not isinstance(series, pd.Series):
        raise ValueError("`series` must be `pd.Series`.")

    plt = ezc.histplot(data=series, bins=bins)
    plt.title = dict(
        text=title,
        subtext=(
            f"Mean: {series.mean().round(round_digits)}. "
            f"Median: {series.median().round(round_digits)}"
        ),
    )
    if x_axis_name is not None:
        plt.xAxis.name = x_axis_name
    if y_axis_name is not None:
        plt.yAxis.name = y_axis_name
    return plt


def tabs_with_histograms(
    data,
    groupby_column,
    data_column,
    plot_title_func,
    x_axis_name,
    y_axis_name,
):
    """Create multiple tabs with a histogram in each.

    Take a `pd.DataFrame`, perform a groupby on the `groupby_column` and then
    generate a tab with a histogram of the values in `data_column` for each group.

    :param data: `pd.DataFrame` containing at least `groupby_column` and `data_column`
    :param groupby_column: column to group by
    :param data_column: column of data to plot in the histogram
    :param plot_title_func: function taking the name of the group and returning the plot
        title
    :param x_axis_name: x axis label
    :param y_axis_name: y axis label
    """
    tabs = Tabs()
    for name, df in data.groupby(groupby_column, observed=True):
        with tabs.add_tab(name):
            title = plot_title_func(name)
            plt = histogram_with_mean_and_median(
                df[data_column].dropna(),
                title=title if title is not None else None,
                x_axis_name=x_axis_name,
                y_axis_name=y_axis_name,
            )
            EZChart(plt, theme=THEME)


def get_summary_table(
    stats_df, flagstat_df, n_reads_total, n_bases_total, secondary=False
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
    n_reads = stats_df.shape[0]
    n_unmapped = flagstat_df["unmapped"].sum()
    n_bases = stats_df["read_length"].sum()
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


def summary(report, sample_names, stats_df, flagstat_df):
    """Create summary section.

    :param report: report object (`ezcharts.components.reports.labs.LabsReport`)
    :param sample_names: collection of sample names
    :param ref_files: collection of reference file names
    :param ref_seqs: collection of reference sequence names
    :param stats_df: `pd.DataFrame` with bamstats per-read stats
    :param flagstat_df: `pd.DataFrame` with bamstats per-file stats
    """
    with report.add_section("Summary", "Summary"):
        # basic stats + intro
        dom_tags.p(
            """
            This report contains visualisations of statistics that can help in
            understanding the results from the wf-human-variation workflow. Each section
            contains different plots or tables, and in general the results are broken
            down by sample. You can quickly jump to an individual section with the links
            in the header bar.
            """
        )
        n_reads_total = stats_df.shape[0]
        n_bases_total = stats_df["read_length"].sum()
        tabs = Tabs()
        for sample_name in [*sample_names]:
            with tabs.add_tab(
                sample_name if sample_name is not None else "total"
            ):
                # get summary stats for individual sample
                get_summary_table(
                    stats_df.query(f"sample_name == '{sample_name}'"),
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


def mapping(report, stats_df_mapped):
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
        with Grid():
            with dom_tags.div():
                sub_heading("Mapping accuracy")
                # Mapping accuracy
                tabs_with_histograms(
                    data=stats_df_mapped,
                    groupby_column="sample_name",
                    data_column="acc",
                    plot_title_func=lambda _: None,
                    x_axis_name="Accuracy [%]",
                    y_axis_name="Number of reads",
                )
            with dom_tags.div():
                sub_heading("Read coverage")
                # coverage per sample
                tabs_with_histograms(
                    data=stats_df_mapped,
                    groupby_column="sample_name",
                    data_column="coverage",
                    plot_title_func=lambda _: None,
                    x_axis_name="Coverage [%]",
                    y_axis_name="Number of reads",
                )


def depths(report, depth_df):
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
            for sample_name, df_samp_file in depth_df.groupby("sample_name"):
                # prepare data for depth vs coordinate plot
                df_depth_vs_coords = (
                    df_samp_file.eval("mean_pos = (start + end) / 2")
                    .eval("step = end - start")
                    .reset_index()
                )
                # Sort column by chromosome number
                df_depth_vs_coords.sort_values(
                    by=["chrom", "start"],
                    key=lambda chrom: chrom.map(CHROMOSOMES)
                    )
                # Reordering
                df_depth_vs_coords['chrom'] = \
                    df_depth_vs_coords.chrom.cat.remove_unused_categories()
                df_depth_vs_coords['chrom'] = \
                    df_depth_vs_coords.chrom.cat.reorder_categories(
                    [i for i in df_depth_vs_coords.chrom.unique()])
                # Extract ref lengths
                ref_lengths = df_depth_vs_coords.groupby(
                    "chrom", observed=True, sort=False)["end"].last()
                total_ref_starts = ref_lengths.cumsum().shift(1, fill_value=0)
                # Add cumulative depth
                df_depth_vs_coords["total_mean_pos"] = df_depth_vs_coords.apply(
                    lambda x: x.mean_pos + total_ref_starts[x.chrom], axis=1)
                # prepare data for cumulative depth plot
                with tabs.add_tab(sample_name):
                    plt = ezc.lineplot(
                        data=df_depth_vs_coords.round(2),
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
