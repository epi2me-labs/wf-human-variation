"""Check (u)BAM files for `@SQ` lines whether they are the same in all headers."""

from pathlib import Path
import sys

import pysam

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    logger = get_named_logger("checkBamHdr")

    if not args.input_path.is_dir():
        raise ValueError(f"Input path '{args.input_path}' must be a directory.")

    target_files = list(args.input_path.glob("*"))
    if not target_files:
        raise ValueError(f"No files found in input directory '{args.input_path}'.")
    # Loop over target files and check if there are `@SQ` lines in all headers or not.
    # Set `is_unaligned` accordingly. If there are mixed headers (either with some files
    # containing `@SQ` lines and some not or with different files containing different
    # `@SQ` lines), set `mixed_headers` to `True`.
    # Also check if there is the SO line, to validate whether the file is (un)sorted.
    first_sq_lines = None
    mixed_headers = False
    sorted_xam = False
    for xam_file in target_files:
        # get the `@SQ` and `@HD` lines in the header
        with pysam.AlignmentFile(xam_file, check_sq=False) as f:
            # compare only the SN/LN/M5 elements of SQ to avoid labelling XAM with
            # same reference but different SQ.UR as mixed_header (see CW-4842)
            sq_lines = [{
                "SN": sq["SN"],
                "LN": sq["LN"],
                "M5": sq.get("M5"),
            } for sq in f.header.get("SQ", [])]
            hd_lines = f.header.get("HD")
        # Check if it is sorted.
        # When there is more than one BAM, merging/sorting
        # will happen regardless of this flag.
        if hd_lines is not None and hd_lines.get('SO') == 'coordinate':
            sorted_xam = True
        if first_sq_lines is None:
            # this is the first file
            first_sq_lines = sq_lines
        else:
            # this is a subsequent file; check with the first `@SQ` lines
            if sq_lines != first_sq_lines:
                mixed_headers = True
                break

    # we set `is_unaligned` to `True` if there were no mixed headers and the last file
    # didn't have `@SQ` lines (as we can then be sure that none of the files did)
    is_unaligned = not mixed_headers and not sq_lines
    # write `is_unaligned` and `mixed_headers` out so that they can be set as env.
    # variables
    sys.stdout.write(
        f"IS_UNALIGNED={int(is_unaligned)};" +
        f"MIXED_HEADERS={int(mixed_headers)};" +
        f"IS_SORTED={int(sorted_xam)}"
    )
    logger.info(f"Checked (u)BAM headers in '{args.input_path}'.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("check_bam_headers_in_dir")
    parser.add_argument("input_path", type=Path, help="Path to target directory")
    return parser
