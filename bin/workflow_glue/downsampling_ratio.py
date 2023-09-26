#!/usr/bin/env python
"""Define if needs downsampling and reports ratio to stdout."""

import sys

from ezcharts.components.mosdepth import load_mosdepth_summary  # noqa: ABS101
import pandas as pd

from .util import wf_parser  # noqa: ABS101


def argparser():
    """Return an arg parser object from arguments."""
    parser = wf_parser("downsampling_ratio")

    parser.add_argument(
        '--summary',
        help="Mosdepth summary file",
        required=True
    )
    parser.add_argument(
        '--bed',
        help="Mosdepth region bed file",
        required=False
    )
    parser.add_argument(
        '--downsample_depth',
        type=float,
        help="Target downsample depth",
        required=True
    )
    parser.add_argument(
        '--margin',
        type=float,
        help="Extra coverage to consider before downsampling",
        required=True
    )

    return parser


def main(args):
    """Run the entry point."""
    # CW-2799: If a region bed is provided, compute the average coverage in the
    # regions only. We leverage the length of each interval by its coverage to reduce
    # the weights of the shorter regions vs the longer.
    if args.bed:
        df = pd.read_csv(args.bed, names=['chrom', 'start', 'end', 'depth'], sep="\t")
        df.eval('length=end-start', inplace=True)
        df.eval('coeff=length*depth', inplace=True)
        value = df.coeff.sum() / df.length.sum()
    # If no bed is provided, simply use the mosdepth summary as a proxy. Use the
    # total_region rather than the total, even though in this case they should be
    # the same.
    else:
        depth_t, depth_r = load_mosdepth_summary(args.summary)
        value = depth_r.loc[depth_r['chrom'] == 'total_region', 'mean'].values[0]
    # Check if the value is at least 10% greater than the expected depth.
    # If it is, compute the downsample ratio. If not, return -1.
    margin = args.margin if args.margin >= 1.0 else 1.0
    if (float(value) / args.downsample_depth) > margin:
        ratio = args.downsample_depth / float(value)
        sys.stdout.write(f"true,{round(ratio, 2)}")
    else:
        sys.stdout.write("false,-1")
