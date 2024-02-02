#!/usr/bin/env python
"""Plot STRs."""

from bokeh.models import BoxZoomTool, ColumnDataSource, HoverTool
from bokeh.models import PanTool, Range1d, ResetTool, WheelZoomTool
from dominate.tags import h3, p, span, table, tbody, td, th, thead, tr
from ezcharts.components import fastcat
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports.labs import LabsReport
from ezcharts.components.theme import LAB_head_resources
from ezcharts.layout.snippets import Grid, Tabs
from ezcharts.plots import BokehPlot
from ezcharts.plots.distribution import histplot
from ezcharts.plots.distribution import MakeRectangles
from natsort import natsorted
import pandas as pd
from .util import wf_parser  # noqa: ABS101


WORKFLOW_NAME = 'wf-human-str'
REPORT_TITLE = f'{WORKFLOW_NAME} report'


def argparser():
    """Return an arg parser object from arguments."""
    parser = wf_parser("report_str")

    parser.add_argument(
        '-o', '--output', required=True, dest="output_report",
        help="Output report")
    parser.add_argument(
        '--params',
        help="workflow parameters"
    )
    parser.add_argument(
        '--versions',
        help="workflow versions"
    )
    parser.add_argument(
        '--vcf',
        help="STRAGLR vcf",
        required=True
    )
    parser.add_argument(
        '--stranger',
        help="STRANGER tsv",
        required=True
    )
    parser.add_argument(
        '--straglr',
        help="STRAGLR tsv",
        required=True
    )
    parser.add_argument(
        '--stranger_annotation',
        help="STRANGER annotated tsv",
        required=True
    )
    parser.add_argument(
        "--sample_name",
        required=True
    )
    parser.add_argument(
        "--str_content",
        required=True,
        help="str content csv"
    )
    parser.add_argument(
        "--read_stats",
        required=True,
        help="read statistics output from bamstats"
    )

    return parser


def add_triangle(plt, x_pos, idx):
    """Draw triangle."""
    plt.add_series(dict(
        type='line',
        datasetIndex=idx,
        markPoint={
            'data': [{
                'symbol': 'triangle',
                'coord': [x_pos, 0],
                'symbolSize': [20, 20],
                'symbolOffset': [0, 15],
                'itemStyle': {
                    'color': 'rgb(255, 255, 255)',
                    'borderColor': 'rgb(0, 0, 0)',
                    'borderWidth': 1.5
                }
            }]
        }
    ))


def add_rectangle(plt, x_start, x_end, height, idx, rectangle):
    """Draw rectangle."""
    if rectangle == 'normal':
        colour = 'rgba(144, 198, 231, 0.4)'
    elif rectangle == 'pathogenic':
        colour = 'rgba(239, 65, 53, 0.4)'
    plt.add_dataset(dict(
        source=[[x_start, x_end, height]],
        dimensions=['x_starts', 'ends', 'heights']
    ))
    plt.add_series(dict(
        type='custom',
        name=rectangle+" range",
        renderItem=MakeRectangles(),
        datasetIndex=idx,
        encode={
            'x': ['x_starts', 'ends'],
            'y': ['heights']
        },
        itemStyle={
            'color': colour},
        clip=True
    ))


def parse_vcf(fname, info_cols=None, nrows=None):
    """Parse a VCF file. The INFO column is parsed to a dictionary.

    :param info_cols: dict of field:dtype for INFO fields to store
        as distinct column.
    :param nrows: number of rows to read from file (including header).
    """
    header = "CHROM POS ID REF ALT QUAL FILTER INFO FORMAT GT".split()
    vcf = pd.read_csv(
        fname, delimiter='\t', comment='#', names=header, nrows=nrows)
    # create a dictionary out of INFO
    vcf['INFO'] = vcf['INFO'].str.split(";") \
        .apply(lambda x: dict([y.split("=") for y in x]))
    # add a column defining the type of the variant
    rlen = vcf['REF'].apply(len)
    alen = vcf['ALT'].apply(len)
    vcf['type'] = 'sub'
    vcf.loc[rlen > alen, 'type'] = 'del'
    vcf.loc[rlen < alen, 'type'] = 'ins'
    # add requested INFO subfields as columns
    if info_cols is not None:
        for field, dtype in info_cols.items():
            try:
                vcf[field] = vcf['INFO'].apply(lambda x: x.get(field, None))
                vcf[field] = vcf[field].astype(dtype)
            except KeyError:
                pass
    return vcf


def process_bamstats(bamstats_file):
    """Load bamstats results into dataframe."""
    d = pd.read_csv(
        bamstats_file, sep="\t", header=0)

    return d


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

    plt = histplot(data=series, bins=bins)
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


def create_str_histogram(
        repeat, hist_data, pathologic_min, pathologic_max, normal_max,
        cn1, cn2, disease):
    """Create a histogram of STR results for a given repeat."""
    h3(disease + ' (' + repeat + ')')

    plt = histplot(data=hist_data[hist_data['VARID'] == repeat]
                   ['copy_number'].values, binwidth=1)
    histogram_data = plt.dataset[0].source
    max_cols = histogram_data.max(axis=0)
    max_rectangle_height = max_cols[2]

    xaxis = {
        'name': "Repeat number",
        'nameGap': '30',
        'nameLocation': 'middle',
        'nameTextStyle': {'fontSize': '14', 'fontStyle': 'bold'},
        'min': '0',
        'max': pathologic_max
    }

    yaxis = {
        'name': "Number of supporting reads",
        'nameGap': '30',
        'nameLocation': 'middle',
        'nameTextStyle': {'fontSize': '14', 'fontStyle': 'bold'},
        'max': max_rectangle_height
    }

    plt.xAxis = xaxis
    plt.yAxis = yaxis

    add_rectangle(
        plt, 0, normal_max, max_rectangle_height, 1, 'normal'
    )
    add_rectangle(
        plt, pathologic_min, pathologic_max, max_rectangle_height, 2,
        'pathogenic'
    )

    # Override plot legend to remove dummy `0` histogram legend entry
    plt.legend = {
        "data": [
            {"name": "normal range"},
            {"name": "pathogenic range"},
        ]}

    add_triangle(plt, cn1, 3)
    add_triangle(plt, cn2, 4)

    EZChart(plt, theme='epi2melabs')


def create_repeat_content_plot(all_data, haplotype_data, repeat_unit):
    """Show repeats and interruptions as rectangular bars."""
    # Creating additional variables for the plot
    unique_reads = haplotype_data['read_id'].unique()
    unique_sequences = all_data['sequence'].unique()
    max_str_sequence = all_data['str_seq_length'].max()
    range_max = max_str_sequence + 1

    # normal_max and pathologic_min are repeat numbers, so get the repeat unit
    # size and use this to convert to length for plotting
    # make a copy of the dataframe to avoid SettingWithCopyWarning
    ru_size = len(repeat_unit)
    new_haplotype_data = haplotype_data.copy()
    new_haplotype_data['str_normal_max'] = new_haplotype_data[
        'str_normal_max'] * ru_size
    new_haplotype_data['str_pathologic_min'] = new_haplotype_data[
        'str_pathologic_min'] * ru_size

    colours = [
        "#1f78b4", "#fdbf6f", "#b2df8a", "#33a02c", "#bbbb88",
        "#baa2a6", "#e08e79", "#B098A4", "#2D4739", "#523249",
        "#AAC0AA", "#120D31", "#F6C0D0", "#F39237", "#CB48B7",
        "#6D9F71", "#337357", "#F49FBC", "#FFD3BA", "#87255B",
        "#B5F44A", "#CCE2A3", "#A5243D", "#F7B538", "#5E747F",
        "#DD1C1A", "#1282A2", "#150578", "#73FBD3", "#D4ADCF",
        "#856084", "#FDE12D", "#FFC0BE", "#B1740F", "#931F1D",
        "#E8C7DE", "#D9BDC5", "#b069a6", "#86467d", "#542c4e",
        "#a2a25d", "#251b65", "#9f2d8e", "#501647", "#ef8fac",
        "#e2366a", "#b31947", "#590d24", "#39f9bf", "#06c68c",
        "#047c58", "#b3ffb3", "#66ff66", "#00e600", "#009900",
        "#004d00", "#ffccff", "#ff4dff", "#e600e6", "#660066"
    ]
    colourmap = {}
    # Mapping colours to unique repeats/interruptions in the seq
    i = 0
    for seq in unique_sequences:
        if seq == all_data['repeat_unit'].any():
            # Consistent colour for the RU's
            colourmap[seq] = "#0173B2"
        else:
            colourmap[seq] = colours[i]
            i += 1
            if i == 60:
                i = 0

    # Setting up data source for plot
    source = ColumnDataSource(
        data=dict(
            read=haplotype_data['read_id'],
            normal_max=haplotype_data['str_normal_max'],
            pathologic_min=haplotype_data['str_pathologic_min'],
            seq=haplotype_data['sequence'],
            truncated_seq=haplotype_data['truncated_seq'],
            seq_length=haplotype_data['length'],
            seq_start=haplotype_data['start'],
            seq_mid=((haplotype_data['end']) - (haplotype_data['length'])/2),
            seq_end=haplotype_data['end'],
            seq_colour=[colourmap[x] for x in haplotype_data['sequence']],
        )
    )
    plt = BokehPlot(
        y_range=list(unique_reads),
        width=1200,
        x_range=Range1d(0, range_max),
        height=40*unique_reads.size,
        outline_line_color=None,  # Stops plot having an outline
    )

    p = plt._fig

    # Add haplotype title
    p.title.text = f"Haplotype: {haplotype_data['haplotype'].values[0]}"

    # Defining which tools to have on the toolbar
    p.toolbar.tools = [
        WheelZoomTool(),
        BoxZoomTool(),
        ResetTool(),
        PanTool()]

    # Adding rectangles to the plot for each repeat/interruption
    p.rect(
        x="seq_mid",  # x-coord of the centre of the rectangle
        y="read",  # y-coord of the centre of the rectangle
        width="seq_length",
        height=0.7,
        source=source,
        fill_alpha=0.6,
        color="seq_colour",
        line_color=None
    )

    # Add normal max line
    p.line(
        x="normal_max",
        y="read",
        line_color="black",
        line_width=1,
        legend_label="Normal Max.",
        source=source
    )

    # Add pathogenic min line
    p.line(
        x="pathologic_min",
        y="read",
        line_color="red",
        line_width=1,
        legend_label="Pathogenic Min.",
        source=source
    )

    p.add_layout(p.legend[0], 'left')

    hover = HoverTool(tooltips=[
        ("Read ID", "@read"),
        ("Index", "@seq_start:@seq_end"),
        ("Sequence", "@truncated_seq")
    ])
    p.add_tools(hover)
    p.yaxis.visible = False

    p.grid.grid_line_color = None

    return plt


def str_content_plot(df, mutation_results):
    """Generate Bokeh plot of read content."""
    for repeat in df['str_identifier'].unique().tolist():
        filtered_df = df[df['str_identifier'] == repeat]
        if (
            filtered_df['varid'].values[0]
            in mutation_results
        ):
            # Adding title and subtext
            h3(repeat)
            repeat_unit_desc = filtered_df['repeat_unit'].values[0]
            p(f"Repeat Unit: {repeat_unit_desc}")

            # Create Bokeh Plot for each haplotype
            haplotypes = filtered_df['haplotype'].unique().tolist()
            haplotypes.sort()

            for haplotype in haplotypes:
                filtered_haplotype_df = filtered_df[
                    filtered_df['haplotype'] == haplotype
                ]
                if not filtered_haplotype_df.empty:
                    # Generate read_count to use for height scaling
                    read_count = len(filtered_haplotype_df['read_id'].unique())

                    # for fewer than 20 reads, pad the plot height to ensure
                    # the plot area is large enough to display them correctly
                    plot_hpx = read_count * 40
                    if read_count < 20:
                        plot_hpx += 50
                    plot_height = f"{plot_hpx}px"

                    # create_report_content_plot() args:
                    # df filtered on STR
                    # df filtered on STR & haplotype
                    EZChart(
                        create_repeat_content_plot(
                            filtered_df,
                            filtered_haplotype_df,
                            repeat_unit_desc
                        ),
                        theme='epi2melabs',
                        height=plot_height)
                p(""" """)  # Adds blank line between plots
        else:
            continue


def make_report(
        read_lengths, params, versions, vcf, stranger, straglr,
        stranger_annotation, str_content_csv, args):
    """Put the STR plots into a report."""
    report = LabsReport(
        f"{args.sample_name} | {REPORT_TITLE}", WORKFLOW_NAME, params,
        versions, head_resources=[*LAB_head_resources])

    vcf_data = parse_vcf(vcf, info_cols={'VARID': str})

    with report.add_section('STR genotyping results', 'STR results'):
        p(
            """
            The tabs below contain short tandem repeat (STR) expansion plots for each
            repeat genotyped in the sample. The coloured boxes denote the normal and
            pathogenic ranges of repeat numbers, and the triangles denote the median
            number of repeats in each allele.
            """
        )

        # set up the results for tabs
        mutation_results = []
        pre_mutation_results = []
        normal_results = []

        # dict for summary table
        summary_data = []

        data = pd.read_csv(straglr, sep="\t", header=0)
        stranger_data = pd.read_csv(stranger, sep="\t", header=0)
        # subtract one from start position to account for how STRs
        # are called by straglr
        stranger_data["POS"] = stranger_data["POS"] - 1
        # read in annotations to get str status
        stranger_annotations = pd.read_csv(stranger_annotation, sep="\t", header=0)
        # subtract one from start position to account for how STRs
        # are called by straglr
        stranger_annotations["POS"] = stranger_annotations["POS"] - 1

        merged = pd.merge(data, stranger_data, left_on='start', right_on='POS')

        for repeat in merged['VARID'].unique().tolist():
            normal_max = merged[merged['VARID'] == repeat][
                'STR_NORMAL_MAX'].unique()
            pathologic_min = merged[merged['VARID'] == repeat][
                'STR_PATHOLOGIC_MIN'].unique()

            str_status = stranger_annotations.loc[stranger_annotations[
                'VARID'] == repeat]['STR_STATUS'].to_string(index=False)
            all_status = str_status.split(",")
            if not (0 < len(all_status) <= 2):
                raise ValueError(f"Invalid number of alleles at site {repeat}")
            # Use the first and last status in the list.
            # If the list report one allele only, it will report the same value
            str_status1 = all_status[0]
            str_status2 = all_status[-1]

            max_copy_number = max(merged[merged['VARID'] == repeat][
                'copy_number'].values)

            disease = merged[merged['VARID'] == repeat]['Disease'].unique()

            repeat_unit = merged[merged['VARID'] == repeat]['repeat_unit'].unique()[0]

            chromosome = merged[merged['VARID'] == repeat]['CHROM'].unique()[0]

            # increment start position by one as this is a variand call and therefore
            # should be 1-based
            start = merged[merged['VARID'] == repeat]['start'].unique()[0] + 1

            end = merged[merged['VARID'] == repeat]['end'].unique()[0]

            # Use the first allele.
            allele1 = merged[merged['VARID'] == repeat]['allele'].unique()[0]

            # Use the last allele (equal to the first if one allele only).
            allele2 = merged[merged['VARID'] == repeat]['allele'].unique()[-1]

            size1 = merged[
                (merged['VARID'] == repeat) &
                (merged['allele'] == allele1)
            ]['size'].median()

            size2 = merged[
                (merged['VARID'] == repeat) &
                (merged['allele'] == allele2)
            ]['size'].median()

            reads = len(merged[merged['VARID'] == repeat])

            if max_copy_number < pathologic_min:
                pathologic_max = pathologic_min + 10
            else:
                pathologic_max = max_copy_number

            genotype = vcf_data.loc[vcf_data['VARID'] == repeat][
                'GT'].to_string(index=False)
            genotype_values = genotype.split(':')
            copy_numbers = genotype_values[2]
            cns = copy_numbers.split('/')
            cn1 = cns[0]
            cn2 = cns[1]

            histogram_vars = [
                repeat, merged, pathologic_min, pathologic_max, normal_max, cn1, cn2,
                disease]

            if 'full_mutation' in str_status:
                mutation_results.append(histogram_vars)
            elif 'pre_mutation' in str_status:
                pre_mutation_results.append(histogram_vars)
            else:
                normal_results.append(histogram_vars)

            a1_badge_colour = ""
            a2_badge_colour = ""

            if str_status1 == "normal":
                a1_badge_colour = "badge bg-primary"
            elif str_status1 == "pre_mutation":
                a1_badge_colour = "badge bg-warning"
            elif str_status1 == "full_mutation":
                a1_badge_colour = "badge bg-danger"

            if str_status2 == "normal":
                a2_badge_colour = "badge bg-primary"
            elif str_status2 == "pre_mutation":
                a2_badge_colour = "badge bg-warning"
            elif str_status2 == "full_mutation":
                a2_badge_colour = "badge bg-danger"

            summary_vars = [
                (repeat, ""),
                (repeat_unit, ""),
                (chromosome, ""),
                (start, ""),
                (end, ""),
                (reads, ""),
                (allele1, a1_badge_colour),
                (size1, ""),
                (allele2, a2_badge_colour),
                (size2, "")
            ]
            summary_data.append(summary_vars)

        sorted_data = natsorted(
            summary_data, key=lambda x: (str(x[2][0]), str(x[3][0])))

        # Set up the tabs
        tabs = Tabs()
        with tabs.add_tab("Mutation", True):
            if not mutation_results:
                p("There are no pathogenic results detected.")
            else:
                for repeat in mutation_results:
                    create_str_histogram(*repeat)
        with tabs.add_tab("Pre-mutation", False):
            if not pre_mutation_results:
                p("There are no pre-mutation results detected.")
            else:
                for repeat in pre_mutation_results:
                    create_str_histogram(*repeat)
        with tabs.add_tab("Normal", False):
            if not normal_results:
                p("There are no normal results detected.")
            else:
                for repeat in normal_results:
                    create_str_histogram(*repeat)

    headers = [
        "Gene", "Repeat unit", "Chromosome", "Start", "End", "No. reads covering locus",
        "Allele 1 repeat number*", "Allele 1 size (bp)*", "Allele 2 repeat number*",
        "Allele 2 size (bp)*"]

    with report.add_section(
            'Summary of STR genotyping results', 'STR summary'):
        p(
            """
            The table below summarises the repeat expansions genotyped in this sample.
            """
        )
        with table(cls="table"):
            with thead():
                for header in headers:
                    th(f"{header}")
            with tbody():
                for record in sorted_data:
                    with tr():
                        for item, badge_colour in record:
                            cell = td()
                            cell.add(span(item, cls=badge_colour))
        p("""* Median values""")
        p("""Repeat number key:""")
        span("Normal", cls="badge bg-primary")
        span("Pre-mutation", cls="badge bg-warning")
        span("Mutation", cls="badge bg-danger")

    with report.add_section('STR content', 'STR content'):
        p(
            """
            The tabs below contain short tandem repeat (STR) expansion plots for each
            repeat genotyped in the sample. These plots display repeat units and
            interruptions for each supporting read. Hover over the coloured bars to view
            more detailed information for each read: the read identifier, position of
            each repeat unit within the read, and the sequence content. Please note,
            repeat content plots are displayed only for those repeats which fall in the
            'Mutation' or 'Pre-mutation' ranges.
            """
        )
        # Load and sort Pandas df
        df_pre_sorting = pd.read_csv(str_content_csv)
        df = df_pre_sorting.sort_values(by=[
            'haplotype',
            'str_seq_length'
        ], ascending=False)
        # Setup results for tabs
        str_mutation_results = []
        str_pre_mutation_results = []
        # Append VARID only to results lists
        for repeat in mutation_results:
            str_mutation_results.append(repeat[0])
        for repeat in pre_mutation_results:
            str_pre_mutation_results.append(repeat[0])

        # Set up the tabs
        tabs = Tabs()
        with tabs.add_tab("Mutation", True):
            if not mutation_results:
                p("There are no pathogenic results detected.")
            else:
                str_content_plot(df, str_mutation_results)

        with tabs.add_tab("Pre-mutation", True):
            if not pre_mutation_results:
                p("There are no pre-mutation results detected.")
            else:
                str_content_plot(df, str_pre_mutation_results)

    with report.add_section('Quality Control Data', 'QC'):
        with Grid():
            max_read_length_to_show = (
                (read_lengths["read_length"] / 1000)
                .quantile(0.99, interpolation="lower")
                .round()
            )

            plt = fastcat.read_length_plot(read_lengths)
            plt.xAxis.max = max_read_length_to_show
            EZChart(plt, theme='epi2melabs')

            plt = histogram_with_mean_and_median(
                read_lengths["mean_quality"],
                title="Mean read quality",
                x_axis_name="Quality",
                y_axis_name="Number of reads"
            )
            EZChart(plt, theme='epi2melabs')

    return report


def main(args):
    """Run the entry point."""
    read_lengths = process_bamstats(args.read_stats)

    report = make_report(
        read_lengths,
        args.params,
        args.versions,
        args.vcf,
        args.stranger,
        args.straglr,
        args.stranger_annotation,
        args.str_content,
        args
    )

    report.write(args.output_report)
