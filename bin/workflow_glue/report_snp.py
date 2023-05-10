#!/usr/bin/env python
"""Create workflow report."""

from aplanat import report
from aplanat.components import bcfstats
from aplanat.components import simple as scomponents
from aplanat.report import WFReport
from .util import wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    report_doc = report.HTMLReport(
        "Haploid variant calling Summary Report",
        ("Results generated through the wf-hap-snp nextflow "
            "workflow provided by Oxford Nanopore Technologies"))

    report_doc = WFReport(
        "Workflow Human SNP", "wf-human-snp",
        revision=args.revision, commit=args.commit)

    section = report_doc.add_section()
    section.markdown("""
### Read Quality control
This section displays basic QC metrics indicating read data quality.
""")
    section = report_doc.add_section()
    bcfstats.full_report(args.vcf_stats, report=section)

    report_doc.add_section(
        section=scomponents.version_table(args.versions))
    report_doc.add_section(
        section=scomponents.params_table(args.params))
    # write report
    report_doc.write(args.report)


def argparser():
    """Create argument parser."""
    parser = wf_parser("report_snp")

    parser.add_argument("report", help="Report output file")
    parser.add_argument(
        "--read_stats", default='unknown',
        help="read statistics output from bamstats")
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--params", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    parser.add_argument(
        "--vcf_stats", default='unknown',
        help="final vcf stats file")

    return parser
