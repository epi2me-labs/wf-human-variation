"""Validate a single (u)BAM file index."""

from pathlib import Path
import sys

import pysam

from .util import get_named_logger, wf_parser  # noqa: ABS101


def validate_xam_index(xam_file):
    """Use fetch to validate the index.

    Invalid indexes will fail the call with a ValueError:
    ValueError: fetch called on bamfile without index
    """
    with pysam.AlignmentFile(xam_file, check_sq=False) as alignments:
        try:
            alignments.fetch()
            has_valid_index = True
        except ValueError:
            has_valid_index = False
    return has_valid_index


def main(args):
    """Run the entry point."""
    logger = get_named_logger("checkBamIdx")

    # Check if a XAM has a valid index
    has_valid_index = validate_xam_index(args.input_xam)
    # write `has_valid_index` out so that they can be set as env.
    sys.stdout.write(
        f"HAS_VALID_INDEX={int(has_valid_index)}"
    )
    logger.info(f"Checked (u)BAM index for: '{args.input_xam}'.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("check_xam_index")
    parser.add_argument("input_xam", type=Path, help="Path to target XAM")
    return parser
