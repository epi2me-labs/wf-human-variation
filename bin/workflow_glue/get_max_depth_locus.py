"""Find max depth window in a `mosdepth` regions BED file and write as locus string."""

from pathlib import Path
import sys

import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    logger = get_named_logger("getMaxDepth")

    # read the regions BED file
    df = pd.read_csv(
        args.depths_bed, sep="\t", header=None, names=["ref", "start", "end", "depth"]
    )

    # get the window with the largest depth
    ref, start, end, depth = df.loc[df["depth"].idxmax()]

    # get the length of the reference of that window
    ref_length = df.query("ref == @ref")["end"].iloc[-1]

    # show the whole reference in case it's shorter than the desired locus size
    if ref_length < args.locus_size:
        start = 1
        end = ref_length
    else:
        # otherwise, show a region of the desired size around the window
        half_size = args.locus_size // 2
        mid = (start + end) // 2
        start = mid - half_size
        end = mid + half_size
        # check if the region starts below `1` or ends beyond the end of the reference
        if start < 1:
            start = 1
            end = args.locus_size
        if end > ref_length:
            start = ref_length - args.locus_size
            end = ref_length

    # write depth and locus string
    sys.stdout.write(f"{depth}\t{ref}:{start}-{end}")

    logger.info("Wrote locus with maximum depth to STDOUT.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("get_max_depth_locus")
    parser.add_argument(
        "depths_bed",
        type=Path,
        help="path to mosdepth regions depth file (can be compressed)",
    )
    parser.add_argument(
        "locus_size", type=int, help="size of the locus in basepairs (e.g. '2000')"
    )
    return parser
