"""Extract unique values for a key from XAM RG DS headers or FASTX comments.

Use pysam to read the description tag of XAM read group header(s) to
collect values for a given key and check the expected cardinality from
one or more XAM; or equivalently the comments of one or more FASTX to
do the same.

Use for example, to extract a single basecaller configuration name in
order to match the input data to suitable models for downstream tools
without troubling the user. No guarantee is made for ordering.
"""

import os
import sys

import pysam

from .util import wf_parser  # noqa: ABS101


# This is not my ideal way to raise nice user-facing errors, but embedding them in the
#  Nextflow process itself is (a) confusing for users; as the error log includes echo
#  commands which may or may not actually be printed and (b) a footgun for developers;
#  who may inadvertently mishandle catching and then (re)returning a bad exit code.
# This also rather neatly keeps intended messaging for users near the code that will
#  cause errors to be raised.
def get_extended_errmsg(key, expected_cardinality):
    """Get an applicable extended error message for a given key and cardinality."""
    if key == "basecall_model" and expected_cardinality == "zero-or-one":
        return """
################################################################################
# INPUT DATA PROBLEM
Your input data contains reads basecalled with more than one basecaller model.

Our workflows automatically select appropriate configuration and models for
downstream tools for a given basecaller model. This cannot be done reliably when
reads with different basecaller models are mixed in the same data set.

## Next steps
To use this workflow you must separate your input files, making sure all reads
are have been basecalled with the same basecaller model.
################################################################################
"""


def path_to_lofn(input_path):
    """Convert the input path to a list of one or more files to be checked."""
    if os.path.isdir(input_path):
        return [
            os.path.join(root, f)
            for (root, dirnames, filenames) in os.walk(input_path)
            for f in filenames
        ]
    else:
        return [input_path]


def xam_extract_ds_key(xam_lofn, key):
    """Extract the set of values for a given key from all RG DS tags."""
    entries = set()
    for xam_fn in xam_lofn:
        with pysam.AlignmentFile(xam_fn, check_sq=False) as xam:
            for read_group in xam.header.get("RG", []):
                for ds_kv in read_group.get("DS", "").split():
                    k, v = ds_kv.split("=", 1)
                    if k == key:
                        entries.add(v)
    return entries


def fastx_extract_ds_key(fastx_lofn, key, stop_after=0):
    """Extract the set of values for a given key from all FASTQ tags."""
    entries = set()
    for fastx_fn in fastx_lofn:
        with pysam.FastxFile(fastx_fn) as fastx:
            for i, read in enumerate(fastx):
                if stop_after > 0 and i > stop_after:
                    break
                for ds_kv in read.comment.split():
                    k, v = ds_kv.split("=", 1)
                    if k == key:
                        entries.add(v)
    return entries


def check_cardinality(obs, desired_cardinality):
    """Return whether the observed cardinality meets the desired cardinality."""
    cardinality_lookup = {
        0: ["zero", "zero-or-one", "zero-or-more"],
        1: ["zero-or-one", "zero-or-more", "one", "one-or-more"],
        2: ["zero-or-more", "one-or-more", "more-than-one"],
    }
    if obs > 1:
        obs = 2
    return desired_cardinality in cardinality_lookup[obs]


def main(args):
    """Script entrypoint.

    Extracts values using the XAM or FASTX extractor and checks the set
    is of the right cardinality. Prints a message to stdout (or stderr)
    and exits appropriately.
    """
    if args.xam:
        extractor = xam_extract_ds_key
        input_path = args.xam
    elif args.fastx:
        extractor = fastx_extract_ds_key
        input_path = args.fastx

    lofn = path_to_lofn(input_path)
    entries = extractor(lofn, args.key)
    if not check_cardinality(len(entries), args.cardinality):
        sys.stdout.write(args.sep.join(entries) + '\n')
        extended_error = get_extended_errmsg(args.key, args.cardinality)
        if args.explode_obviously and extended_error:
            sys.stderr.write(extended_error)
        else:
            sys.stderr.write(
                f"Required {args.cardinality} {args.key} but found {len(entries)}\n"
            )
        sys.exit(os.EX_DATAERR)
    sys.stdout.write(args.sep.join(entries) + '\n')
    sys.exit(os.EX_OK)


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("get_ds_records")
    parser.add_argument("--key", required=True)
    input_arg = parser.add_mutually_exclusive_group(required=True)
    input_arg.add_argument("--xam", help="Path to a single XAM or folder of XAM")
    input_arg.add_argument("--fastx", help="Path to a single FASTX or folder of FASTX")
    parser.add_argument(
        "--sep",
        default="\n",
        help=(
            "Value separator to use if more than one element is printed to stdout."
        ),
    )
    parser.add_argument(
        "--cardinality",
        default="zero-or-more",
        choices=[
            "zero",
            "zero-or-one",
            "zero-or-more",
            "one",
            "one-or-more",
            "more-than-one",
        ],
        help=(
            "Expected cardinality of entries. Will exit EX_DATAERR if the wrong "
            "number of elements are found. Defaults to zero-or-more."
        ),
    )
    parser.add_argument(
        "--explode_obviously",
        action="store_true",
        help=(
            "If appropriate, print a more fulsome user facing error if the expected"
            "cardinality has been violated."
        )
    )
    return parser
