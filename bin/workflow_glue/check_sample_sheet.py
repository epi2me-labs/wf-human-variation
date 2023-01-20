"""Script to check that sample sheet is well-formatted."""
import sys

import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run entry point."""
    logger = get_named_logger("check-sheet")

    try:
        logger.info(f"Reading {args.sample_sheet}.")
        samples = pd.read_csv(args.sample_sheet, sep=None)
        if 'alias' in samples.columns:
            if 'sample_id' in samples.columns:
                sys.stderr.write(
                    "Warning: sample sheet contains both 'alias' and "
                    'sample_id, using the former.')
            samples['sample_id'] = samples['alias']
        if not set(['sample_id', 'barcode']).intersection(samples.columns):
            raise IOError()
    except Exception:
        raise IOError(
            "Could not parse sample sheet, it must contain two columns "
            "named 'barcode' and 'sample_id' or 'alias'.")
    # check duplicates
    dup_bc = samples['barcode'].duplicated()
    dup_sample = samples['sample_id'].duplicated()
    if any(dup_bc) or any(dup_sample):
        raise IOError(
            "Sample sheet contains duplicate values.")
    samples.to_csv(args.output, sep=",", index=False)
    logger.info(f"Written cleaned-up sheet to {args.output}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("check-sample-sheet")
    parser.add_argument('sample_sheet')
    parser.add_argument('output')
    return parser
