"""Check if a sample sheet is valid."""
import csv
import sys

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    logger = get_named_logger("checkSheet")

    barcodes = []
    aliases = []
    sample_types = []
    allowed_sample_types = [
        "test_sample", "positive_control", "negative_control", "no_template_control"
        ]

    try:
        with open(args.sample_sheet, "r") as f:
            csv_reader = csv.DictReader(f)
            n_row = 0
            for row in csv_reader:
                n_row += 1
                if n_row == 1:
                    n_cols = len(row)
                else:
                    # check we got the same number of fields
                    if len(row) != n_cols:
                        raise ValueError(
                            f"Unexpected number of cells in row number {n_row}."
                        )
                try:
                    barcodes.append(row["barcode"])
                except KeyError:
                    sys.stdout.write("'barcode' column missing")
                    sys.exit()
                try:
                    aliases.append(row["alias"])
                except KeyError:
                    sys.stdout.write("'alias' column missing")
                    sys.exit()
                try:
                    sample_types.append(row["type"])
                except KeyError:
                    pass
    except Exception as e:
        sys.stdout.write(f"Parsing error: {e}")
        sys.exit()

    # check barcode and alias values are unique
    if len(barcodes) > len(set(barcodes)):
        sys.stdout.write("values in 'barcode' column not unique")
        sys.exit()
    if len(aliases) > len(set(aliases)):
        sys.stdout.write("values in 'alias' column not unique")
        sys.exit()

    if sample_types:
        # check if "type" column has unexpected values
        unexp_type_vals = set(sample_types) - set(allowed_sample_types)

        if unexp_type_vals:
            sys.stdout.write(
                f"found unexpected values in 'type' column: {unexp_type_vals}. "
                f"Allowed values are: {allowed_sample_types}"
            )
            sys.exit()

        if args.required_sample_types:
            for required_type in args.required_sample_types:
                if required_type not in allowed_sample_types:
                    sys.stdout.write(f"Not an allowed sample type: {required_type}")
                    sys.exit()
                if sample_types.count(required_type) < 1:
                    sys.stdout.write(
                        f"Sample sheet requires at least 1 of {required_type}")
                    sys.exit()

    logger.info(f"Checked sample sheet {args.sample_sheet}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("check_sample_sheet")
    parser.add_argument("sample_sheet", help="Sample sheet to check")
    parser.add_argument(
        "--required_sample_types",
        help="List of required sample types. Each sample type provided must "
             "appear at least once in the sample sheet",
        nargs="*"
    )
    return parser
