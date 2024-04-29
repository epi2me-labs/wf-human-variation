#!/usr/bin/env python
"""Create SNV report."""

import json
import os
import sys

from dominate.tags import a, h6, p
from ezcharts.components.bcfstats import load_bcfstats
from ezcharts.components.clinvar import load_clinvar_vcf
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports.labs import LabsReport
from ezcharts.components.theme import LAB_head_resources
from ezcharts.layout.snippets import DataTable, Stats, Tabs
from ezcharts.plots import util
from ezcharts.plots.categorical import barplot
from ezcharts.plots.matrix import heatmap
import numpy as np
import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101

# Global variables
Colors = util.Colors

# Mutation profile palette to match the COSMIC
# patterns.
cmap = {
    'C>A': Colors.cerulean,
    'C>G': Colors.black,
    'C>T': Colors.cinnabar,
    'T>A': Colors.grey70,
    'T>C': Colors.medium_spring_bud,
    'T>G': Colors.fandango,
}


def indel_sizes(df):
    """Extract indel sizes with added 0-values."""
    df['nlength'] = df['length (deletions negative)'].astype(int)
    df['count'] = df['number of sites'].astype(int)
    # pad just to pull out axes by a minimum
    counts = df.groupby('nlength') \
        .agg(count=pd.NamedAgg(column='count', aggfunc='sum')) \
        .reset_index()
    # Create min/max range of values
    all_counts = pd.DataFrame({
        'nlength': np.arange(counts.nlength.min() - 1, counts.nlength.max() + 1)
    })
    all_counts = pd.merge(all_counts, counts, on="nlength", how="left").fillna(0)
    # Add values where needed:
    return all_counts


def parse_changes(bcftools_dt):
    """Parse changes from counts using same method as in aplanat."""
    df = bcftools_dt
    sim_sub = {
        'G>A': 'C>T', 'G>C': 'C>G', 'G>T': 'C>A',
        'T>A': 'A>T', 'T>C': 'A>G', 'T>G': 'A>C'}

    def canon_sub(sub):
        b1 = sub[0]
        if b1 not in {'A', 'C'}:
            return canon_sub(sim_sub[sub])
        else:
            return b1, sub[2]

    df['canon_sub'] = df['type'].apply(canon_sub)
    df['Reference allele'] = df['canon_sub'].apply(lambda x: x[0])
    df['Alternative allele'] = df['canon_sub'].apply(lambda x: x[1])
    df['count'] = df['count'].astype(int)
    df = df[['Reference allele', 'Alternative allele', 'count']] \
        .groupby(['Reference allele', 'Alternative allele']) \
        .agg(count=pd.NamedAgg(column='count', aggfunc='sum')) \
        .reset_index()
    df2 = df.pivot(
        index="Alternative allele",
        columns="Reference allele",
        values="count")
    return df2


def main(args):
    """Run the entry point."""
    logger = get_named_logger("report_snp")
    clinvar_docs_url = "https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/"

    # Check that the input files exist
    if not os.path.exists(args.vcf_stats):
        raise FileNotFoundError(f"File {args.vcf_stats} not found.")

    # Load the data
    try:
        bcfstats = load_bcfstats(
            args.vcf_stats,
            sample_names=[args.sample_name])
    except IndexError:
        bcfstats = {'SN': pd.DataFrame(), 'TSTV': pd.DataFrame()}

    # Save values to an output JSON
    with open(f'{args.sample_name}.snvs.json', 'w') as json_file:
        json.dump(
            {
                'SNVs': int(bcfstats['SN']['SNPs'].values[0]),
                'Indels': int(bcfstats['SN']['indels'].values[0]),
                'Transition/Transversion rate': bcfstats['TSTV']['ts/tv'].values[0]
            }, json_file)

    if args.skip_report:
        logger.info('Skip report')
        sys.exit(0)
    # Instantiate the report
    report = LabsReport(
        "Small variation statistics", "wf-human-variation",
        args.params, args.versions,
        head_resources=[*LAB_head_resources])

    # VCF At-a-glance report

    with report.add_section('At a glance', 'Summary'):
        tabs = Tabs()
        with tabs.add_tab(args.sample_name):
            if bcfstats['SN'].empty:
                p('No Summary Numbers (SN) lines to plot from the bcftools stats file.')
            else:
                bcfstats['SN'].columns = [
                    i.replace('SNP', 'SNV').replace('MNP', 'MNV')
                    for i in bcfstats['SN'].columns
                ]
                titv = bcfstats['TSTV']['ts/tv'].values[0]
                nsites = bcfstats['SN']['records'].values[0]
                nsnvs = bcfstats['SN']['SNVs'].values[0]
                nindels = bcfstats['SN']['indels'].values[0]
                Stats(
                    columns=4,
                    items=[
                        (f'{"{:,}".format(int(nsites))}', 'Variants'),
                        (f'{"{:,}".format(int(nsnvs))}', 'SNVs'),
                        (f'{"{:,}".format(int(nindels))}', 'Indels'),
                        (f'{titv}', 'Ti/Tv')
                    ])

    # Base statistics
    with report.add_section('Statistics', 'Stats'):
        tabs = Tabs()
        with tabs.add_tab(args.sample_name):
            if bcfstats['SN'].empty:
                p('No Summary Numbers (SN) lines to plot from the bcftools stats file.')
            else:
                DataTable.from_pandas(
                    bcfstats['SN'].drop(columns=['id']),
                    use_index=False)
                DataTable.from_pandas(
                    bcfstats['TSTV'].drop(columns='id'),
                    use_index=False)

    # ClinVar variants
    if not args.skip_annotation:
        if args.clinvar_vcf is not None:
            if os.path.exists(args.clinvar_vcf):
                with report.add_section('ClinVar variant annotations', 'ClinVar'):
                    p(
                        "The ",
                        a("SnpEff", href="https://pcingola.github.io/SnpEff/"),
                        " annotation tool has been used to annotate with",
                        a("ClinVar", href="https://www.ncbi.nlm.nih.gov/clinvar/"), '.'
                        " Variants with ClinVar annotations will appear in the ",
                        "table below, ranked according to their significance. ",
                        "'Pathogenic', 'Likely pathogenic', and 'Unknown ",
                        "significance' will be displayed first, in that order. ",
                        "Please note that variants classified as 'Benign' or ",
                        "'Likely benign' are not reported in this table, but ",
                        "will appear in the VCF output by the workflow. For further",
                        " details on the terms in the 'Significance' column, please",
                        " visit ",
                        a("this page", href=clinvar_docs_url),
                        '.')
                    # check if there are any ClinVar sites to report
                    clinvar_for_report = load_clinvar_vcf(args.clinvar_vcf)
                    if clinvar_for_report.empty:
                        h6('No ClinVar sites to report.')
                    else:
                        DataTable.from_pandas(
                            clinvar_for_report, export=True, use_index=False)

    else:
        # Annotations were skipped
        with report.add_section('ClinVar variant annotations', 'ClinVar'):
            p(
                "This report was generated without annotations. To see"
                " them, re-run the workflow without --skip_annotation.")

    # Change type
    with report.add_section('Substitution types', 'Types'):
        p('Base substitutions aggregated across all samples (symmetrised by pairing).')
        tabs = Tabs()
        with tabs.add_tab(args.sample_name):
            if bcfstats['ST'].empty:
                p(
                    'No Substitution Types (ST) lines to plot from the bcftools '
                    'stats file.'
                )
            else:
                subtype = parse_changes(bcfstats['ST'])
                plt = heatmap(subtype)
                EZChart(plt, 'epi2melabs')

    # Plot mutation spectra if provided
    with report.add_section('Indels length', 'Indels'):
        p('Insertion and deletion lengths aggregated across all samples.')
        tabs = Tabs()
        with tabs.add_tab(args.sample_name):
            if bcfstats['IDD'].empty:
                p(
                    'No Indels Distribution (IDD) lines to plot from the bcftools '
                    'stats file.'
                )
            else:
                sizes = indel_sizes(bcfstats['IDD'])
                plt = barplot(data=sizes, x="nlength", y="count", color=Colors.cerulean)
                EZChart(plt, 'epi2melabs')

    # write report
    report.write(args.report)
    logger.info(f"Written report to '{args.report}'.")


def argparser():
    """Create argument parser."""
    parser = wf_parser("report")

    parser.add_argument("report", help="Report output file")
    parser.add_argument(
        "--read_stats", default='unknown',
        help="Read statistics output from bamstats"
    )
    parser.add_argument(
        "--vcf_stats", default='unknown',
        help="Final VCF stats file"
    )
    parser.add_argument(
        "--clinvar_vcf", required=True,
        help="VCF file of variants annotated in ClinVar"
    )
    parser.add_argument(
        "--sample_name", default='Sample',
        help="Sample name"
    )
    parser.add_argument(
        "--skip_report",
        action='store_true',
        required=False,
        help="Skip generation of HTML report and only output JSON metrics"
    )
    parser.add_argument(
        "--versions", required=True,
        help="Directory containing CSVs containing name,version"
    )
    parser.add_argument(
        "--params", required=True,
        help="Directory containing workflow parameters"
    )
    parser.add_argument(
        "--skip_annotation", action="store_true",
        help="Do not show ClinVar variants in report"
    )

    return parser
