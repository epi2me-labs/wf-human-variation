#!/usr/bin/env python
"""Create workflow report."""

import json

from dominate.tags import a, h4, h6, p
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports.labs import LabsReport
from ezcharts.components.theme import LAB_head_resources
from ezcharts.layout.snippets import DataTable, Grid, Stats, Tabs
from ezcharts.layout.snippets import Progress
from ezcharts.layout.util import isolate_context
from ezcharts.plots import util
from ezcharts.plots.categorical import barplot
from ezcharts.plots.ideogram import ideogram
import numpy as np
import pandas as pd

from .report_utils import read_data  # noqa: ABS101
from .report_utils.common import CHROMOSOMES  # noqa: ABS101
from .report_utils.common import CLINVAR_BASE, NCBI_GENE_BASE  # noqa: ABS101
from .util import wf_parser  # noqa: ABS101

vcf_cols = [
    'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
    'INFO', 'FORMAT', 'TUMOR', 'NORMAL'
    ]
Colors = util.Colors


def read_vcf(fname):
    """Read input VCF as pandas dataframe."""
    new_types = {
        'CHROM': str,
        'SVLEN': int,
        'END': int,
    }
    try:
        df = pd.read_csv(fname, comment='#', header=None, delimiter='\t')
        df = df.rename(columns=dict(enumerate(vcf_cols)))
        df = expand_column(df)
        df = df.drop('RNAMES', 1)
        for (key, dtype) in new_types.items():
            df[key] = df[key].astype(dtype)
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
    return vcf_df.groupby('SVTYPE').agg(**{
        'Count': ('POS', 'count'),
        'Min. Length': ('SVLEN', lambda x: np.min(x.abs())),
        'Ave. Length': ('SVLEN', lambda x: np.median(x.abs())),
        'Max. Length': ('SVLEN', lambda x: np.max(x.abs()))}).transpose()


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
                Stats(
                    columns=3,
                    items=[
                        (f'{"{:,}".format(inserts.shape[0])}',
                         'Number of Insertions'),
                        (f'{"{:,}".format(delets.shape[0])}',
                         'Number of Deletions'),
                        (f'{"{:,}".format(vcf_df.CHROM.unique().shape[0])}',
                         'Chromosomes with SVs'),
                    ])


def sv_size_plots(vcf_data):
    """Plot size distributions of SV calls per type."""
    tabs = Tabs()
    for (index, sample_name, vcf_df) in vcf_data:
        with tabs.add_tab(sample_name):
            if vcf_df.empty:
                p("The workflow found no structural variants to report.")
            else:
                p("This section shows the size distributions of SV calls per type. "
                    "Deletions have negative values.")
                with Grid():
                    # Extract insertions
                    inserts = vcf_df.loc[vcf_df['SVTYPE'] == 'INS']
                    if not inserts.empty:
                        inserts = inserts[['SVTYPE', 'SVLEN']] \
                            .groupby('SVLEN') \
                            .count() \
                            .reset_index()
                        inserts.columns = ['Length', 'Count']
                        # Define max value
                        maxins = int(inserts.Length.max()) \
                            if inserts.Length.max() else 0
                        # Create min/max range of values
                        inserts = pd.merge(pd.DataFrame({
                            'Length': range(0, maxins + 10)
                        }), inserts, on="Length", how="left").fillna(0)
                        # Add insertion plot
                        plt = barplot(data=inserts, x="Length", y="Count")
                        plt.title = {"text": "Insertion size distribution"}
                        EZChart(plt, 'epi2melabs')
                    else:
                        h4("No insertions to plot")
                    # Extract deletions
                    delets = vcf_df.loc[vcf_df['SVTYPE'] == 'DEL']
                    if not delets.empty:
                        delets = delets[['SVTYPE', 'SVLEN']] \
                            .groupby('SVLEN') \
                            .count() \
                            .reset_index()
                        delets.columns = ['Length', 'Count']
                        # Define max value
                        mindel = int(delets.Length.min()) \
                            if delets.Length.min() else 0
                        # Create min/max range of values
                        delets = pd.merge(pd.DataFrame({
                            'Length': range(mindel - 10, 1)
                        }), delets, on="Length", how="left").fillna(0)
                        # Add deletion plot
                        plt = barplot(data=delets, x="Length", y="Count")
                        plt.title = {"text": "Deletion size distribution"}
                        EZChart(plt, 'epi2melabs')
                    else:
                        h4("No deletions to plot")


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
                df = vcf_df[['CHROM', 'POS', 'ID', 'SVLEN']]
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

    # Input all ClinVar VCFs
    clinvar_vcf_data = []
    for index, sample_vcf in enumerate(args.clinvar_vcf):
        # Annotation is on
        if not args.skip_annotation:
            vcf_df = read_data.parse_vcf(sample_vcf)
            vcf_is_empty = read_data.parse_vcf_for_size(sample_vcf)
            clinvar_vcf_data.append((
                index, sample_vcf.split('.')[0], vcf_data, vcf_is_empty[0]))
        else:
            continue

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

    if not args.skip_annotation:
        with report.add_section('ClinVar variant annotations', 'ClinVar'):
            clinvar_docs_url = "https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/"
            p(
                "The ",
                a("SnpEff", href="https://pcingola.github.io/SnpEff/"),
                " annotation tool has been used to annotate with",
                a("ClinVar", href="https://www.ncbi.nlm.nih.gov/clinvar/"), '.'
                " If any variants have ClinVar annotations, they will appear in a ",
                "table below. Please note, this table excludes variants with 'Benign'",
                " or 'Likely benign' significance, however these variants will appear",
                "in the VCF output by the workflow. For further details on the terms",
                " in the 'Significance' column, please visit ",
                a("this page", href=clinvar_docs_url),
                '.')
            tabs = Tabs()
            for (index, sample_name, vcf_data, vcf_is_empty) in clinvar_vcf_data:
                with tabs.add_tab(sample_name):
                    if vcf_is_empty:
                        h6('No ClinVar sites were found.')
                    else:
                        # display ClinVar records
                        data_table = DataTable(
                            headers=[
                                'Chrom',
                                'Pos',
                                'Gene(s)',
                                'ClinVar',
                                'Significance',
                                'Type',
                                'Consequence',
                                'HGVSg'],
                            export=True,
                            file_name=(
                                f'{sample_name}-wf-human-variation-sv-clinvar'))

                        for variant in vcf_data:
                            significance = ", ".join(variant.info['CLNSIG'])
                            benigns = ["Benign", "benign"]
                            if any([x in significance for x in benigns]):
                                continue
                            else:
                                try:
                                    all_ncbi_urls = []
                                    clinvar_gene_string = variant.info['GENEINFO']
                                    all_genes = clinvar_gene_string.split('|')
                                    for gene in all_genes:
                                        gene_symbol, gene_id = gene.split(':')
                                        ncbi_url = f"""
                                            <a href=
                                                \"{NCBI_GENE_BASE}{gene_id}\">{
                                                    gene_symbol}</a>"""
                                        all_ncbi_urls.append(ncbi_url)
                                    ncbi_gene_url = ", ".join(all_ncbi_urls)
                                except KeyError:
                                    ncbi_gene_url = "No affected genes found"

                                clinvar_id = variant.id
                                with isolate_context():
                                    clinvar_url = a(
                                        clinvar_id, href='%s%s' % (
                                            CLINVAR_BASE, clinvar_id))

                                significance = significance.replace(
                                    "_", " ").capitalize()

                                variant_type = variant.info['CLNVC']
                                variant_type = variant_type.replace(
                                    "_", " ").capitalize()
                                if variant_type == 'Single nucleotide variant':
                                    variant_type = 'SNV'

                                consequences = []
                                try:
                                    all_consequences = variant.info['MC']
                                    for each_conseq in all_consequences:
                                        ontology, consequence = each_conseq.split('|')
                                        consequence = consequence.replace(
                                            "_", " ").capitalize()
                                        consequence = consequence.replace(
                                            "5 prime utr", "5' UTR")
                                        consequence = consequence.replace(
                                            "3 prime utr", "3' UTR")
                                        consequences.append(consequence)
                                    consequences = ", ".join(consequences)
                                except KeyError:
                                    consequences = "No consequences found"

                                hgvs = variant.info['CLNHGVS']
                                hgvs = ", ".join(hgvs)

                                data_table.add_row(
                                    title=None,
                                    columns=[
                                        variant.chrom,
                                        variant.pos,
                                        ncbi_gene_url,
                                        clinvar_url,
                                        significance,
                                        variant_type,
                                        consequences,
                                        hgvs])
    else:
        # Annotations were skipped
        with report.add_section('ClinVar variant annotations', 'ClinVar'):
            p(
                "This report was generated without annotations. To see"
                " them, re-run the workflow without --skip_annotation.")

    with report.add_section('Karyogram', 'Karyogram'):
        karyoplot(vcf_data, args)

    with report.add_section('Size distribution', 'Size'):
        sv_size_plots(vcf_data)

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
        help="Report output file.")
    parser.add_argument(
        "--vcf",
        nargs='+',
        required=True)
    parser.add_argument(
        "--clinvar_vcf",
        nargs='+',
        required=False)
    parser.add_argument(
        "--genome",
        default='hg38',
        required=False)
    parser.add_argument(
        "--eval_results",
        nargs='+',
        required=False)
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    parser.add_argument(
        "--params", default=None, required=True,
        help="A csv containing the parameter key/values")
    parser.add_argument(
        "--params-hidden", default="",
        help="Comma delimited list of keys to hide from parameters table")
    parser.add_argument(
        "--versions", required=True,
        help="directory contained CSVs containing name,version.")
    parser.add_argument(
        "--skip_annotation", action="store_true",
        help="Do not show ClinVar variants in report.")

    return parser
