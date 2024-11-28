#!/usr/bin/env python
"""check_seq_ref.

Compare a FASTA reference to BAM/CRAM SQ lines and determine
if realignment is required.
"""

import os
import sys

import pysam
from .util import wf_parser  # noqa: ABS101

# NOTE Both OK and DATAERR are permissible exits from this script so
#      if you want to raise an error to stop Nextflow, it better be
#      something else!


def argparser():
    """Create argument parser."""
    parser = wf_parser("check_sq_ref")

    parser.add_argument(
        "--xam",
        required=True,
        help="Input BAM or CRAM"
    )
    parser.add_argument(
        "--ref",
        required=True,
        help="Reference FASTA file"
    )
    return parser


def main(args):
    """Run entry point."""
    try:
        xam = pysam.AlignmentFile(args.xam, check_sq=False)
        ref = pysam.FastaFile(args.ref)
    except ValueError:
        sys.stderr.write(
            "[FAIL] One (or both) of the input files could not be"
            " read. Are they the right format?"
        )
        sys.exit(os.EX_NOINPUT)

    ref_reflen = set(zip(ref.references, ref.lengths))
    xam_reflen = set(zip(xam.references, xam.lengths))
    diff = ref_reflen ^ xam_reflen

    if len(diff) > 0:
        sys.stdout.write(" ".join([
            "sequence_name",
            "sequence_length",
            "in_ref",
            "in_xam",
        ]) + '\n')
        for reflen in diff:
            sys.stdout.write(" ".join([
                f"{reflen[0]}",
                f"{reflen[1]}",
                '1' if reflen in ref_reflen else '0',
                '1' if reflen in xam_reflen else '0',
            ]) + '\n')

        # honestly flake8 why
        sys.stderr.write(
            "[FAIL] There is at least one (name, length)"
            " pair that does not map 1:1 between the input alignment"
            " and reference.\n"
        )
        sys.exit(os.EX_DATAERR)
    else:
        sys.stderr.write(
            "[OKAY] Input alignment and reference sequences map 1:1\n"
        )
