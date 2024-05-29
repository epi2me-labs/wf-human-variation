#!/usr/bin/env python
"""Check whether the input is a modbam."""

import os
import sys

import pysam

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    # Load other input files
    logger = get_named_logger("check_valid_modbam")
    logger.info(f'Checking file: {args.bam}')

    # Check the first 10K reads of the input bam for ML/MM fields
    valid_reads = 0
    fields = ['mm', 'ml']
    for i, alignment in enumerate(pysam.AlignmentFile(args.bam)):
        n_tags = len([
            tag for (tag, val) in alignment.get_tags() if tag.lower() in fields])
        if n_tags == 2:
            valid_reads += 1
            break
        if i >= 9999:
            break
    if valid_reads == 0:
        sys.exit(os.EX_DATAERR)


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("check_modbam")
    parser.add_argument("bam", help="Input bam file")
    return parser
