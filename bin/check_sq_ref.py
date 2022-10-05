#!/usr/bin/env python
"""check_seq_ref.

Compare a FASTA reference to BAM/CRAM SQ lines and determine
if realignment is required.
"""

import os
import sys

import pysam

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--xam", "--bam", "--cram",
        required=True,
    )
    parser.add_argument(
        "--ref",
        required=True,
    )
    args = parser.parse_args()

    xam = pysam.AlignmentFile(args.xam)
    ref = pysam.FastaFile(args.ref)

    ref_reflen = set(zip(ref.references, ref.lengths))
    xam_reflen = set(zip(xam.references, xam.lengths))
    diff = ref_reflen ^ xam_reflen

    if len(diff) > 0:
        for reflen in diff:
            print(" ".join([
                "sequence_name",
                "sequence_length",
                "in_ref",
                "in_xam",
            ]))
            print(
                reflen[0],
                reflen[1],
                '1' if reflen in ref_reflen else '0',
                '1' if reflen in xam_reflen else '0',
            )

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
