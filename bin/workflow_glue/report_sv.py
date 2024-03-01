#!/usr/bin/env python
"""Create workflow report."""

import json

from dominate.tags import a, h4, p
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports.labs import LabsReport
from ezcharts.components.theme import LAB_head_resources
from ezcharts.layout.snippets import DataTable, Stats, Tabs
from ezcharts.layout.snippets import Progress
from ezcharts.plots import util
from ezcharts.plots.distribution import histplot
from ezcharts.plots.ideogram import ideogram
import numpy as np
import pandas as pd


from .report_utils.common import CHROMOSOMES  # noqa: ABS101
from .util import wf_parser  # noqa: ABS101

vcf_cols = [
    'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
    'INFO', 'FORMAT', 'TUMOR', 'NORMAL'
    ]
Colors = util.Colors


def read_vcf(fname):
    """Read input VCF as pandas dataframe."""
    # CW-2492: BND have no END/SVLEN flags, so we do conversion to int/flow
    # contextually when needed
    try:
        df = pd.read_csv(fname, comment='#', header=None, delimiter='\t')
        df = df.rename(columns=dict(enumerate(vcf_cols)))
        df = expand_column(df)
        df = df.drop(columns='RNAMES').astype({'CHROM': str})
        df['CHROM'] = df['CHROM'].str.replace('chr', '')
        df['CHROM'] = df['CHROM'].astype(str)
        df['CHROM'] = df['CHROM'].str.replace('chr', '')
        df = df.loc[df['CHROM'].isin(CHROMOSOMES)]
    except pd.errors.EmptyDataError:
        df = pd.DataFrame()
    return df


def expand_column(df, column='INFO', row_split_delim=';', key_val_delim='='):
    """Expand a column for a given dataframe."""
    data = []
    for row in df[column].items():
        split_row = row[1].split(row_split_delim)
        split_row_items = {}
        for i in split_row:
            key_val_split = i.split(key_val_delim)
            # Replace SVINSLEN with SVLEN since referring to the same metric
            if key_val_split[0] == 'SVINSLEN':
                key_val_split[0] = 'SVLEN'
            if not len(key_val_split) > 1:
                split_row_items[key_val_split[0]] = ''
                continue
            split_row_items[key_val_split[0]] = key_val_split[1]

        data.append(split_row_items)

    df2 = pd.DataFrame(data)
    df2.fillna('', inplace=True)
    return pd.merge(df, df2, how='left', left_index=True, right_index=True)


def get_sv_summary_table(vcf_df):
    """Aggregate summary info for SV calls per type."""
    # CW-2492: replace missing values with 0s to allow stats, and ensure
    # int type for SVTYPE
    return vcf_df\
        .replace('', 0)\
        .astype({'SVLEN': int})\
        .groupby('SVTYPE')\
        .agg(**{
            'Count': ('POS', 'count'),
            'Min. Length': ('SVLEN', lambda x: np.min(x.abs())),
            'Ave. Length': ('SVLEN', lambda x: np.median(x.abs())),
            'Max. Length': ('SVLEN', lambda x: np.max(x.abs()))})\
        .transpose()


def sv_stats(vcf_data):
    """Display base stats."""
    tabs = Tabs()
    for (index, sample_name, vcf_df) in vcf_data:
        with tabs.add_tab(sample_name):
            if vcf_df.empty:
                p("The workflow found no structural variants to report.")
            else:
                inserts = vcf_df.loc[vcf_df['SVTYPE'] == 'INS']
                delets = vcf_df.loc[vcf_df['SVTYPE'] == 'DEL']
                # CW-2492: count other SV types as well
                other = vcf_df\
                    .loc[(vcf_df['SVTYPE'] != 'INS') & (vcf_df['SVTYPE'] != 'DEL')]
                Stats(
                    columns=4,
                    items=[
                        (f'{"{:,}".format(inserts.shape[0])}',
                         'Number of Insertions'),
                        (f'{"{:,}".format(delets.shape[0])}',
                         'Number of Deletions'),
                        (f'{"{:,}".format(other.shape[0])}',
                         'Other SV types'),
                        (f'{"{:,}".format(vcf_df.CHROM.unique().shape[0])}',
                         'Chromosomes with SVs'),
                    ])
                p('Other SV types are: inversions, duplications and translocations.')


def centered_bins(dataset, values, binning_mode='sqrt'):
    """Create bins distributed around the zero value."""
    # Compute bins
    edges = np.histogram_bin_edges(
        dataset[values].dropna(), bins=binning_mode, range=None, weights=None)
    bin_sizes = [edges[n+1]-edges[n] for n in range(0, len(edges)-1)]
    bin_size = int(np.ceil(max(bin_sizes)))

    # Define min and max values from the dataset
    min_val = int(dataset[values].values.min())
    max_val = int(dataset[values].values.max())

    # If min_val >= 0, no deletions to consider
    if min_val >= 0:
        del_itvs = []
    else:
        del_itvs = [
            -i for i in range(bin_size, abs(min_val) + bin_size, bin_size)][::-1]

    # If max_val <= 0, no insertions to consider
    if max_val <= 0:
        ins_itvs = []
    else:
        ins_itvs = list(range(bin_size, max_val + bin_size, bin_size))

    # Return all bins
    return del_itvs + [0] + ins_itvs


def sv_size_plots(vcf_data, max_size=5000):
    """Plot size distributions of SV calls per type."""
    tabs = Tabs()
    for (index, sample_name, vcf_df) in vcf_data:
        with tabs.add_tab(sample_name):
            if vcf_df.empty:
                p("The workflow found no structural variants to report.")
            else:
                p("This section shows the size distributions of SV calls per type. "
                    "Deletions have negative values.")
                # Extract deletions
                # CW-2492: ensure SVLEN is int
                indels = vcf_df\
                    .loc[(vcf_df['SVTYPE'] == 'DEL') | (vcf_df['SVTYPE'] == 'INS')]\
                    .astype({'SVLEN': int, 'END': int})
                # Keep SVs within range of interest
                indels = indels.loc[np.abs(indels['SVLEN'].values) < max_size]
                # Create plot
                if not indels.empty:
                    # Define bin intervals
                    bins = centered_bins(indels, 'SVLEN')
                    indels = indels['SVLEN']
                    indels.columns = ['Length']
                    # Add deletion plot
                    plt = histplot(data=indels, bins=bins, stat='count')
                    # override excharts axisLabel interval
                    plt.xAxis = dict(
                        name='Length',
                        axisLabel=dict(
                            interval="auto",
                            rotate=30
                        ),
                        max=max_size,
                        min=-max_size
                        )
                    plt.title = {"text": "Indels size distribution"}
                    EZChart(plt, 'epi2melabs')
                    p("The plot shows Indels with |length| < 5Kb.")
                else:
                    h4("No Indels to plot")


def karyoplot(vcf_data, args):
    """Karyogram plot."""
    p("Chromosomal hotspots of structural variation.")
    tabs = Tabs()
    for (index, sample_name, vcf_df) in vcf_data:
        with tabs.add_tab(sample_name):
            if vcf_df.empty:
                p("The workflow found no structural variants to report.")
            else:
                colors = {
                    "INS": Colors.cinnabar,
                    "DEL": Colors.cerulean
                }
                # CW-2492: ensure SVLEN is int
                df = vcf_df\
                    .loc[(vcf_df['SVTYPE'] == 'INS') | (vcf_df['SVTYPE'] == 'DEL')]\
                    .astype({'SVLEN': int})
                df = df[['CHROM', 'POS', 'ID', 'SVLEN']]
                df.columns = ['chr', 'start', 'name', 'length']
                # Define end point to start+length...
                df = df.eval('end=start + length')
                # ... but if its too small it won't appear, so set it
                # to at least 100Kb in size to be able to show the hotspots
                df.loc[df['length'] < 1e5, 'end'] += 1e5 - \
                    df.loc[df['length'] < 1e5, 'length']
                # Define color mappings
                df['color'] = vcf_df['SVTYPE'].map(colors).fillna(Colors.black)
                df = df[['chr', 'start', 'end', 'name', 'color']]
                # Prepare the ideogram
                plt = ideogram(blocks=df, genome=args.genome)
                EZChart(plt, height='600px', width='90%', theme='epi2melabs')
                p("""Red: Insertion""")
                p("""Blue: Deletion""")


def main(args):
    """Run the entry point."""
    # Input all VCFs
    vcf_data = []
    for index, sample_vcf in enumerate(args.vcf):
        vcf_df = read_vcf(sample_vcf)
        vcf_data.append((index, sample_vcf.split('.')[0], vcf_df))

    # Create report file
    report = LabsReport(
        "Structural variants analysis", "wf-human-variation",
        args.params, args.versions,
        head_resources=[*LAB_head_resources])

    with report.add_section('At a glance', 'Summary'):
        p(
            "This section displays a description"
            ' of the variant calls made by ',
            a("Sniffles", href="https://github.com/fritzsedlazeck/Sniffles"), '.')
        sv_stats(vcf_data)

    with report.add_section('Variant calling results', 'Variants'):
        p(
            "This section displays summary statistics"
            " of the variant calls made by",
            a("Sniffles", href="https://github.com/fritzsedlazeck/Sniffles"), '.')
        tabs = Tabs()
        for (index, sample_name, vcf_df) in vcf_data:
            with tabs.add_tab(sample_name):
                if vcf_df.empty:
                    p("The workflow found no structural variants to report.")
                else:
                    DataTable.from_pandas(get_sv_summary_table(vcf_df))

    with report.add_section('Karyogram', 'Karyogram'):
        karyoplot(vcf_data, args)

    with report.add_section('Size distribution', 'Size'):
        sv_size_plots(vcf_data, args.hist_sv_max_length)

    # Import jsons and create the progress bars before plotting
    values = []
    if args.eval_results:
        for (index, sample_name, vcf_df) in vcf_data:
            with open(args.eval_results[index]) as f:
                data = json.load(f)
                tp = data['TP-call']
                fp = data['FP']
                fn = data['FN']
                # Represent the metrics as progress bars ranging from 0 to 1
                f1 = Progress(
                    value_min=0.0,
                    value_max=1.0,
                    value_now=round(data['f1'], 2),
                    bar_cls="progress-bar-striped",
                    height=50)
                pr = Progress(
                    value_min=0.0,
                    value_max=1.0,
                    value_now=round(data['precision'], 2),
                    bar_cls="progress-bar-striped",
                    height=50)
                rc = Progress(
                    value_min=0.0,
                    value_max=1.0,
                    value_now=round(data['recall'], 2),
                    bar_cls="progress-bar-striped",
                    height=50)
                # Append values
                values.append({
                    'TP': tp, 'FP': fp, 'FN': fn,
                    'F1': f1, 'PR': pr, 'RC': rc
                    })

    with report.add_section('Results evaluation', 'Evaluation'):
        if not args.eval_results:
            p(
                "This report was generated without evaluation"
                " results. To see them, re-run the workflow with"
                " --sv_benchmark set.")
        else:
            p(
                "This section displays the ",
                a("Truvari", href="https://github.com/ACEnglish/truvari"),
                " benchmarking and"
                " evaluation of the variant calls made by",
                a("Sniffles", href="https://github.com/fritzsedlazeck/Sniffles"), '.')
            tabs = Tabs()
            for (index, sample_name, vcf_df) in vcf_data:
                with tabs.add_tab(sample_name):
                    #
                    # Evaluation results
                    #
                    vals = values[index]
                    p(
                        "This sections displays the truvari"
                        " evaluation metrics for your SV calls.")
                    # Show statistics as 3x2 matrix of statistics.
                    Stats(
                        columns=3,
                        items=[
                            (vals['TP'], 'True positives'),
                            (vals['FP'], 'False positives'),
                            (vals['FN'], 'False negatives'),
                            (vals['F1'], 'F1-score'),
                            (vals['PR'], 'Precision'),
                            (vals['RC'], 'Recall')
                        ])

    #
    # write report
    #
    report.write(args.output)


def argparser():
    """Create argument parser."""
    parser = wf_parser("report_sv")
    parser.add_argument(
        "output",
        help="Report output file."
    )
    parser.add_argument(
        "--vcf",
        nargs='+',
        required=True
    )
    parser.add_argument(
        "--genome",
        default='hg38',
        required=False
    )
    parser.add_argument(
        "--hist_sv_max_length",
        default=5000,
        help="Max length of an SV to display in the histogram",
        required=False
    )
    parser.add_argument(
        "--eval_results",
        nargs='+',
        required=False
    )
    parser.add_argument(
        "--revision", default='unknown',
        help="Git branch/tag of the executed workflow"
    )
    parser.add_argument(
        "--commit", default='unknown',
        help="Git commit of the executed workflow"
    )
    parser.add_argument(
        "--params", default=None, required=True,
        help="A CSV containing the parameter key/values"
    )
    parser.add_argument(
        "--params-hidden", default="",
        help="Comma delimited list of keys to hide from parameters table"
    )
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version."
    )

    return parser
