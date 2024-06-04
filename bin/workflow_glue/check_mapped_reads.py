#!/usr/bin/env python
"""check_mapped_reads.

Check that the first read of a BAM file is mapped.
"""

import os
import sys

import pysam
from .util import wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("check_mapped_reads")

    parser.add_argument(
        "--xam",
        required=True,
        help="Input BAM or CRAM"
    )
    return parser


def main(args):
    """Run entry point."""
    try:
        xam = pysam.AlignmentFile(args.xam, check_sq=False)
    except ValueError:
        sys.stderr.write(
            "[FAIL] The input file could not be read. Is it in the right format?"
        )
        sys.exit(os.EX_NOINPUT)

    for aln in xam:
        sys.stdout.write(f"has_maps={int(aln.is_mapped)};")
        sys.exit()
