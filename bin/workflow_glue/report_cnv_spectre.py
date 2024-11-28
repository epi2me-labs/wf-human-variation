#!/usr/bin/env python
"""Plot Spectre CNVs."""

from dominate.tags import a, p
from ezcharts.components.reports.labs import LabsReport
from ezcharts.components.theme import LAB_head_resources
from ezcharts.layout.snippets import DataTable
from ezcharts.layout.snippets import Stats
from ezcharts.plots import util
import pandas as pd

from .report_utils.common import CHROMOSOMES  # noqa: ABS101
from .util import wf_parser  # noqa: ABS101

# Setup simple globals
WORKFLOW_NAME = 'wf-human-variation'
REPORT_TITLE = f'{WORKFLOW_NAME} CNV report'
Colors = util.Colors


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report_cnv_spectre")
    parser.add_argument(
        '--sample_id',
        help="sample_id"
    )
    parser.add_argument(
        '--params',
        help="workflow parameters"
    )
    parser.add_argument(
        '--versions',
        help="workflow versions"
    )
    parser.add_argument(
        '--cnv_bed',
        help="BED file of CNV calls"
    )
    parser.add_argument(
        "--workflow_version", required=True,
        help="Workflow version",
    )
    parser.add_argument(
        "--karyotype",
        help="Text file containing karyotype"
    )
    parser.add_argument(
        '-o', '--output', required=True, dest="output_report",
        help="Output report")

    return parser


def make_report(params, versions, cnv_df, args):
    """Create the report."""
    report = LabsReport(
        f"{args.sample_id} | {REPORT_TITLE}",
        WORKFLOW_NAME, params, versions, args.workflow_version,
        head_resources=[*LAB_head_resources])

    with open(args.karyotype) as karyotype_file:
        line = karyotype_file.readline().rstrip()
        if line:
            karyotype = line
        else:
            karyotype = "Not applicable"

    if not cnv_df.empty:
        # get the totals
        total_cnv = cnv_df.shape[0]
        type_counts = cnv_df['type'].value_counts()
        total_dups = type_counts.get('DUP', 0)
        total_dels = type_counts.get('DEL', 0)

        spectre_url = "https://github.com/nanoporetech/ont-spectre"

        with report.add_section('At a glance', 'Summary'):
            p(
                "This section displays a summary of the CNVs detected using an ",
                a("ONT implementation of Spectre", href=spectre_url),
                ", as part of the wf-human-variation workflow."
            )
            Stats(
                columns=4,
                items=[
                    (karyotype, 'Predicted karyotype*'),
                    (total_cnv, 'Total no. of CNVs'),
                    (total_dups, 'Total no. of duplications'),
                    (total_dels, 'Total no. of deletions')
                ])

            p("""*The predicted karyotype is "Not applicable" if chromosomes X and
                Y are absent from the input data or have insufficient coverage.
                """)

        with report.add_section(
                'Spectre CNV calls', 'Spectre CNV calls'):

            p("""This table lists the CNVs called by Spectre. A link to Ensembl
            genes within each CNV region is provided. However, this is not available for
            copy number events greater than 5Mb, as this is the size limit for the query
            used to generate the link.""")

            data_table = DataTable(
                headers=[
                    'Chr',
                    'Start',
                    'End',
                    'Type',
                    'Size',
                    'Ensembl genes'])

            for index, row in cnv_df.iterrows():
                if row['size'] > 5000000:
                    genes_in_region = "Unavailable"
                else:
                    ensembl = "https://rest.ensembl.org/overlap/region/human/"
                    region = f"{row['chr']}:{row['start']}-{row['end']}"
                    options = "?feature=gene;content-type=text/x-gff3"
                    genes_in_region = f"""
                        <a href=\"{ensembl}{region}{options}\">Genes</a>"""
                data_table.add_row(
                    title=None,
                    columns=[
                        row['chr'],
                        row['start'],
                        row['end'],
                        row['type'],
                        row['size'],
                        genes_in_region])
    else:
        with report.add_section(
                'Spectre CNV calls', 'Spectre CNV calls'):
            p("""No CNVs were detected.""")

    return report


def main(args):
    """Run the entry point."""
    # read CNV BED into df
    try:
        cnv_df = pd.read_csv(args.cnv_bed, delim_whitespace=True, header=None)
        cnv_df.columns = ['chr', 'start', 'end', 'type', 'size', 'call']

        # hg19 chromosomes are read in as integers, so convert to string to ensure
        # filtering is successful
        filtered_df = cnv_df[cnv_df['chr'].astype(str).isin(CHROMOSOMES)]
    except pd.errors.EmptyDataError:
        filtered_df = pd.DataFrame()

    report = make_report(
        args.params,
        args.versions,
        filtered_df,
        args)

    report.write(args.output_report)
