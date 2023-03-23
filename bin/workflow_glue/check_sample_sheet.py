"""Check if a sample sheet is valid."""
import csv
import sys

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    logger = get_named_logger("checkSheet")

    barcodes = []
    aliases = []
    types = []

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
                    exit()
                try:
                    aliases.append(row["alias"])
                except KeyError:
                    sys.stdout.write("'alias' column missing")
                    exit()
                try:
                    types.append(row["type"])
                except KeyError:
                    pass
    except Exception as e:
        sys.stdout.write(f"Parsing error: {e}")
        exit()

    # check barcode and alias values are unique
    if len(barcodes) > len(set(barcodes)):
        sys.stdout.write("values in 'barcode' column not unique")
        exit()
    if len(aliases) > len(set(aliases)):
        sys.stdout.write("values in 'alias' column not unique")
        exit()

    if types:
        # check if "type" column has unexpected values
        unexp_type_vals = set(types) - set(
            [
                "test_sample",
                "positive_control",
                "negative_control",
                "no_template_control",
            ]
        )
        if unexp_type_vals:
            sys.stdout.write(
                f"found unexpected values in 'type' column: {unexp_type_vals}. "
                "allowed values are: `['test_sample', 'positive_control', "
                "'negative_control', 'no_template_control']`"
            )
            exit()

    logger.info(f"Checked sample sheet {args.sample_sheet}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("check_sample_sheet")
    parser.add_argument("sample_sheet", help="Sample sheet to check")
    return parser
