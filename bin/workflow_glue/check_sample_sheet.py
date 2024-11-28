"""Check if a sample sheet is valid."""
import codecs
import csv
import os
import re
import sys

from .util import get_named_logger, wf_parser  # noqa: ABS101


# Some Excel users save their CSV as UTF-8 (and occasionally for a reason beyond my
# comprehension, UTF-16); Excel then adds a byte order mark (unnecessarily for UTF-8
# I should add). If we do not handle this with the correct encoding, the mark will
# appear in the parsed data, causing the header to be malformed.
# See CW-2310
def determine_codec(f):
    """Peek at a file and return an appropriate reading codec."""
    with open(f, 'rb') as f_bytes:
        # Could use chardet here if we need to expand codec support
        initial_bytes = f_bytes.read(8)

        for codec, encoding_name in [
            [codecs.BOM_UTF8, "utf-8-sig"],  # use the -sig codec to drop the mark
            [codecs.BOM_UTF16_BE, "utf-16"],  # don't specify LE or BE to drop mark
            [codecs.BOM_UTF16_LE, "utf-16"],
            [codecs.BOM_UTF32_BE, "utf-32"],  # handle 32 for completeness
            [codecs.BOM_UTF32_LE, "utf-32"],  # again skip LE or BE to drop mark
        ]:
            if initial_bytes.startswith(codec):
                return encoding_name
        return None  # will cause file to be opened with default encoding


def main(args):
    """Run the entry point."""
    logger = get_named_logger("checkSheet")

    barcodes = []
    aliases = []
    sample_types = []
    analysis_groups = []
    allowed_sample_types = [
        "test_sample", "positive_control", "negative_control", "no_template_control"
    ]

    if not os.path.exists(args.sample_sheet) or not os.path.isfile(args.sample_sheet):
        sys.stdout.write("Could not open sample sheet file.")
        sys.exit()

    try:
        encoding = determine_codec(args.sample_sheet)
        with open(args.sample_sheet, "r", encoding=encoding) as f:
            try:
                # Excel files don't throw any error until here
                csv.Sniffer().sniff(f.readline())
                f.seek(0)  # return to initial position again
            except Exception as e:
                # Excel fails with UniCode error
                sys.stdout.write(
                    "The sample sheet doesn't seem to be a CSV file.\n"
                    "The sample sheet has to be a CSV file.\n"
                    "Please verify that the sample sheet is a CSV file.\n"
                    f"Parsing error: {e}"
                 )

                sys.exit()

            csv_reader = csv.DictReader(f)
            n_row = 0
            for row in csv_reader:
                n_row += 1
                if n_row == 1:
                    n_cols = len(row)
                else:
                    # check we got the same number of fields
                    if len(row) != n_cols:
                        sys.stdout.write(
                            f"Unexpected number of cells in row number {n_row}"
                        )
                        sys.exit()
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
                try:
                    analysis_groups.append(row["analysis_group"])
                except KeyError:
                    pass
    except Exception as e:
        sys.stdout.write(f"Parsing error: {e}")
        sys.exit()

    # check barcodes are correct format
    for barcode in barcodes:
        if not re.match(r'^barcode\d\d+$', barcode):
            sys.stdout.write("values in 'barcode' column are incorrect format")
            sys.exit()

    # check barcodes are all the same length
    first_length = len(barcodes[0])
    for barcode in barcodes[1:]:
        if len(barcode) != first_length:
            sys.stdout.write("values in 'barcode' column are different lengths")
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
    if analysis_groups:
        # if there was a "analysis_group" column, make sure it had values for all
        # samples
        if not all(analysis_groups):
            sys.stdout.write(
                "if an 'analysis_group' column exists, it needs values in each row"
            )
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
