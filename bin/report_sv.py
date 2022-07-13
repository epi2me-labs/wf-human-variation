#!/usr/bin/env python
"""Create workflow report."""

import argparse
from datetime import datetime
import json
import sys

import aplanat
from aplanat import bio, hist, report
from aplanat.components import depthcoverage, fastcat
import aplanat.graphics
from bokeh.layouts import gridplot, layout
from bokeh.models import BasicTickFormatter, Range1d
import numpy as np
import pandas as pd


vcf_cols = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']


def read_vcf(fname):
    """Read a VCF file."""
    df = pd.read_csv(fname, comment='#', header=None, delimiter='\t')
    df = df.rename(columns=dict(enumerate(vcf_cols)))
    return df


def expand_column(df, column='INFO', row_split_delim=';', key_val_delim='='):
    """Expand a column for a given dataframe."""
    data = []
    for row in df[column].items():
        split_row = row[1].split(row_split_delim)
        split_row_items = {}
        for i in split_row:
            key_val_split = i.split('=')
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


def get_sv_karyograms(vcf_df, sv_types, sv_colours):
    """Plot positions of SV calls per type."""
    karyograms = list()
    for sv, col in zip(sv_types, sv_colours):
        data = vcf_df.loc[vcf_df['SVTYPE'] == sv[0]]
        if data.empty:
            continue
        plot = bio.karyotype(
            [data['POS']],
            [data['CHROM']],
            [sv[1]],
            [col],
            alpha=0.2,
            height=300)
        aplanat.export_jsx(plot, f'{sv[0]}_karyogram.jsx')
        karyograms.append(plot)
    return karyograms


def get_sv_size_plots(vcf_df, sv_types, sv_colours):
    """Plot size distributions of SV calls per type."""
    length_plots = list()
    for sv, col in zip(sv_types, sv_colours):
        data = np.log10(vcf_df.loc[vcf_df['SVTYPE'] == sv[0], 'SVLEN'].abs())
        if data.empty:
            continue
        plot = hist.histogram(
            [data], colors=[col],
            names=[sv[1]], bins=200, xlim=(1, None),
            height=250,
            title="{} SV lengths".format(sv[1]),
            x_axis_label='log10(SV length / bases)',
            y_axis_label='count')
        aplanat.export_jsx(plot, f'{sv[0]}_size.jsx')
        length_plots.append(plot)
    return length_plots


def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser(sys.argv[1:])
    parser.add_argument(
        "output",
        help="Report output file.")
    parser.add_argument(
        "--vcf",
        nargs='+',
        required=True)
    parser.add_argument(
        "--reads_summary",
        nargs='+',
        required=False)
    parser.add_argument(
        "--eval_results",
        nargs='+',
        required=False)
    parser.add_argument(
        "--read_depth",
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
    args = parser.parse_args()
    report_doc = report.WFReport(
        "wf-human-sv report", "wf-human-sv",
        revision=args.revision, commit=args.commit)

    for index, sample_vcf in enumerate(args.vcf):
        sample_name = sample_vcf.split('.')[0]

        #
        # Front matter
        #
        section = report_doc.add_section()
        section.markdown(
            f"## Sample: {sample_name}")
        section.markdown(
            f"```Date: {datetime.today().strftime('%Y-%m-%d')}```")

        #
        # Input dataset QC
        #
        sample_summary = None
        for summary in args.reads_summary or []:
            if summary.split('.')[0] == sample_name:
                sample_summary = summary

        for read_depth in args.read_depth:
            if read_depth.split('.')[0] == sample_name:
                read_depth = read_depth

        if sample_summary:
            reads_summary_df = pd.read_csv(sample_summary, sep='\t')
            read_qual = fastcat.read_quality_plot(reads_summary_df)
            read_length = fastcat.read_length_plot(reads_summary_df)
            read_length.x_range = Range1d(0, 100000)
            read_length.xaxis.formatter = BasicTickFormatter(
                use_scientific=False)
            section = report_doc.add_section()
            section.markdown("### Read Quality Control")
            section.markdown(
                "This sections displays basic QC"
                " metrics indicating read data quality.")
            section.plot(layout(
                [[read_length, read_qual]],
                sizing_mode="stretch_width")
            )
            rd_plot = depthcoverage.cumulative_depth_from_dist(read_depth)
            section = report_doc.add_section()
            section.markdown("### Genome coverage")
            section.markdown(
                "Genome wide read coverage")
            section.plot(gridplot(
                [rd_plot],
                ncols=2))

        #
        # Variant calls
        #
        section = report_doc.add_section()
        section.markdown("### Variant calling results")
        section.markdown(
            "This section displays a summary view"
            " of the variant calls made by cuteSV.")

        # Assuming we are using hg37 or 38 here
        chroms_37 = [str(x) for x in range(1, 23)] + ['X', 'Y']
        sv_types = (('INS', 'Insertion'), ('DEL', 'Deletion'))
        sv_colours = ['red', 'green']

        vcf_df = read_vcf(sample_vcf)
        vcf_df = expand_column(vcf_df)
        vcf_df = vcf_df.drop('RNAMES', 1)
        vcf_df['CHROM'] = vcf_df['CHROM'].astype(str)
        vcf_df['SVLEN'] = vcf_df['SVLEN'].astype(int)
        vcf_df['CHROM'] = vcf_df['CHROM'].str.replace('chr', '')
        vcf_df = vcf_df.loc[vcf_df['CHROM'].isin(chroms_37)]

        table = get_sv_summary_table(vcf_df)
        karyograms = gridplot(
            get_sv_karyograms(vcf_df, sv_types, sv_colours),
            ncols=2)
        size_plots = gridplot(
            get_sv_size_plots(vcf_df, sv_types, sv_colours),
            ncols=2)

        section.table(
            table,
            index=True,
            sortable=False,
            paging=False,
            searchable=False)
        section.plot(
            layout(
                [karyograms],
                [size_plots],
                sizing_mode="stretch_width"))

        if not args.eval_results:
            section.markdown(
                "This report was generated without evaluation"
                " results. To see them, re-run the workflow with"
                " --mode benchmark set.")
            continue

        #
        # Evaluation results
        #
        section.markdown("### Evaluation results")
        data = None
        with open(args.eval_results[index]) as f:
            data = json.load(f)
        section = report_doc.add_section()
        section.markdown(
            "This sections displays the truvari"
            " evaluation metrics for your SV calls.")
        exec_summary = aplanat.graphics.InfoGraphItems()
        exec_summary.append(
            "TP-call",
            "{:.2f}".format(data['TP-call']),
            "chart-pie"
        )
        exec_summary.append(
            "FP",
            "{:.2f}".format(data['FP']),
            "chart-pie"
        )
        exec_summary.append(
            "FN",
            "{:.2f}".format(data['FN']),
            "chart-pie"
        )
        exec_summary.append(
            "F1",
            "{:.2f}".format(data['f1']),
            "chart-pie"
        )
        exec_summary.append(
            "Precision",
            "{:.2f}".format(data['precision']),
            "chart-pie"
        )
        exec_summary.append(
            "Recall",
            "{:.2f}".format(data['recall']),
            "chart-pie"
        )
        exec_plot = aplanat.graphics.infographic(
            exec_summary.values(), ncols=4)
        section.plot(exec_plot, key="exec-plot")

    #
    # Params reporting
    #
    section = report_doc.add_section()
    section.markdown("## Workflow parameters")
    section.markdown(
        "The table below highlights values of"
        " the main parameters used in this analysis.")
    params = []
    hidden_params = args.params_hidden.split(',')
    with open(args.params) as f:
        params_data = json.load(f)
    for key, value in params_data.items():
        if key not in hidden_params:
            params.append((key, value))
    df_params = pd.DataFrame(params, columns=['Key', 'Value'])
    section.table(
        df_params, sortable=False, paging=False, index=False, searchable=False)

    #
    # Software versions
    #
    section = report_doc.add_section()
    section.markdown("## Software versions")
    section.markdown('''The table below highlights versions
                    of key software used within the analysis''')
    versions = list()
    if args.versions is not None:
        with open(args.versions) as fh:
            for line in fh.readlines():
                name, version = line.strip().split(',')
                versions.append((name, version))
    versions = pd.DataFrame(versions, columns=('Name', 'Version'))
    section.table(versions, index=False)

    #
    # write report
    #
    report_doc.write(args.output)


if __name__ == "__main__":
    main()
