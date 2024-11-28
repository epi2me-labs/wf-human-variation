#!/usr/bin/env python
"""Unify VCFs for tertiary analyses."""

import gzip
import os

from .util import get_named_logger, wf_parser  # noqa: ABS101


# Get logger
logger = get_named_logger("report_sv")


def get_lines(file_path):
    """
    Return lines in a file.

    This is done in order to handle both zipped and unzipped files.
    """
    try:
        with gzip.open(file_path, "rt") as file_h:
            lines = file_h.readlines()
    except gzip.BadGzipFile:
        with open(file_path) as file_h:
            lines = file_h.readlines()
    return lines


def modify_repeat_line(repeat_line):
    """
    Add SVTYPE=REP to INFO in the vcf.

    This is necessary for repeat vcf for providers that are not ONT, to enable
    Geneyx application reading these lines as REP effect.
    ONT repeat files already contain SVTYPE=STR in their INFO section, which is
    also recognized by Geneyx as REP effect
    """
    (chrom, pos, vid, ref, alt, qual, vfilter, info, vformat, _) = (
        repeat_line.split("\t")
    )
    info = info + ";SVTYPE=REP"
    repeat_line = "\t".join(
        (chrom, pos, vid, ref, alt, qual, vfilter, info, vformat, _)
    )
    return repeat_line


def is_valid_line(line):
    """
    Check that the line contains SVTYPE=.

    this is Geneyx way to recognizing structural variants
    also checks that the gt (genotype) is not "./." (which means a no-call)
    a line is valid if it contains SVTYPE= and it had a valid genotype
    """
    # check if line is header or null / empty
    if line.startswith("#") or not line:
        return True

    cells = line.split("\t")
    info_cell = cells[7]
    format_key = cells[8]
    format_value = cells[9]

    # if the "svtype" line is exists the line valid
    if "SVTYPE=" in info_cell:
        return True

    gt_index = format_key.split(":").index("GT")
    gt_value = format_value.split(":")[gt_index]

    # if "svtype" not present and GT is 0 or "./." the line is invalid
    if gt_index == "0" or gt_value == "./.":
        return False

    return True


def write_vcf_content(ftype, lines_dict, output_h, skip_svtype):
    """
    Write VCF.

    Writes the content of a file into the output file skips the header
    of the file modifies lines when necessary
    """
    start_read_content_flag = False
    for vcf_line in lines_dict[ftype]:
        if start_read_content_flag:
            if ftype == "repeat" and not skip_svtype:
                vcf_line = modify_repeat_line(vcf_line)
            output_h.write(vcf_line)
        if vcf_line.startswith("#CHROM"):
            start_read_content_flag = True

    if not start_read_content_flag:
        logger.warning(f"No vcf header for given {ftype} file")


def get_vcf_header_from_lines(vcf_lines):
    """Get VCF lines and returns the header."""
    vcf_header = ""
    for line in vcf_lines:
        vcf_header = vcf_header + line
        if line.find("#CHROM") != -1:
            break

    return vcf_header


def is_empty(lines):
    """Check if the file is empty."""
    for line in lines:
        striped_line = line.strip()
        if striped_line != "" and not striped_line.startswith("#"):
            return False

    return True


def combine_headers(files_lines):
    """Combine the headers of multiple VCF files."""
    c_line = None
    h_lines = []
    printed_format = False

    # Process each file header
    for key in ['sv', 'cnv', 'repeat']:
        lines = files_lines[key]
        if lines:
            for line in lines:
                if line.startswith('##') and line not in h_lines:
                    if line.startswith('##fileformat') and not printed_format:
                        h_lines.append(line)
                        printed_format = True
                        continue
                    elif not line.startswith('##fileformat'):
                        h_lines.append(line)
                elif line.startswith('#CHROM') and not c_line:
                    c_line = line
    h_lines.append(c_line)
    return h_lines


def create_unified_file(files_lines, output_path, skip_svtype):
    """
    Create the unified file.

    Prints the sv file complete (with header), then the cnv without a header
    and repeats without a header and with modified info section
    if sv file is empty prints all cnv lines (with header) and repeat
    lines without a header if cnv is also empty prints repeat lines with the
    header
    """
    # We need to combine the headers first, to avoid issues downstream when
    # sorting/indexing woth bcftools and keep the output compliant.
    header = combine_headers(files_lines)

    # Then, combine the body of the files
    with open(output_path, "w+") as output_h:
        # Save the combined header first
        output_h.write("".join(header))

        # prints into the output file all the sv file
        if files_lines["sv"] is not None:
            write_vcf_content("sv", files_lines, output_h, skip_svtype)

        # adds the CNV lines (without the header)
        if files_lines["cnv"] is not None:
            write_vcf_content("cnv", files_lines, output_h, skip_svtype)

        # Then add the repeat file
        if files_lines["repeat"] is not None:
            write_vcf_content("repeat", files_lines, output_h, skip_svtype)

        output_h.flush()
        output_h.close()


# creates a dictionary with the different vcf files per type and
# calls the unifying function
def main(args):
    """Run the entry point."""
    # Define input variables.
    output_path = args.outputPath
    sv_path = args.svPath
    cnv_path = args.cnvPath
    repeat_path = args.repeatPath
    skip_svtype = True
    output_file_path = output_path

    # Check that there is at least one input
    if not sv_path and not cnv_path and not repeat_path:
        logger.warning("Empty/no vcf files to concatenate")
        return

    struct_files = {
        "sv": sv_path,
        "cnv": cnv_path,
        "repeat": repeat_path
    }
    struct_lines = {}

    for file_type in struct_files.keys():
        struct_lines[file_type] = None

        if struct_files[file_type] is None:
            logger.warning(
                f'File type "{file_type}" not provided.'
            )
            continue

        # This is not in the original UnifyVcf script, but better check if
        # the file is there.
        if not os.path.isfile(struct_files[file_type]):
            logger.warning(
                f'File "{struct_files[file_type]}" does not exists.'
            )
            continue

        lines = get_lines(struct_files[file_type])
        if is_empty(lines):
            logger.warning(
                f'Can\'t unify file type "{file_type}". '
                "The file is empty."
            )
            continue

        struct_lines[file_type] = lines

    create_unified_file(struct_lines, output_file_path, skip_svtype)


def argparser():
    """Create argument parser."""
    parser = wf_parser("UnifyVcf")
    parser.add_argument(
        '-o', '--outputPath',
        help='the unified output VCF path (required)',
        required=True
    )
    parser.add_argument(
        '-s', '--svPath',
        help='SV input file path (optional)',
        required=False, default=None
    )
    parser.add_argument(
        '-c', '--cnvPath',
        help='CNV input file path (optional)',
        required=False, default=None
    )
    parser.add_argument(
        '-r', '--repeatPath',
        help='repeats input file path (optional)',
        required=False, default=None
    )
    return parser
