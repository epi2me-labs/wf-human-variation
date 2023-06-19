#!/usr/bin/env python
"""Fix malformed QDNAseq VCFs by reformatting files with a single CNV call.

See CW-1491 for more information.
"""

import argparse
import sys


def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser("fix_vcf")
    parser.add_argument(
        '-i', required=True, dest="vcf",
        help="Input VCF")
    parser.add_argument(
        '-o', required=True, dest="fixed_vcf",
        help="Fixed VCF")
    parser.add_argument(
        '--sample_id', required=True, dest="sample_id")
    return parser


def main(args):
    """Run entry point."""
    result = open(args.fixed_vcf, 'w')

    seen_header = False
    in_record = False
    current_record = []
    with open(args.vcf, 'r') as vcf_file:
        for line_i, line in enumerate(vcf_file):
            # pass metaheaders through
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    seen_header = True
                result.write(line)
            else:
                if not seen_header:
                    # reached a data line without seeing a header, write it out
                    sys.stderr.write(
                        "[fix_vcf] Missing header, will insert one to fix.\n")
                    result.write(
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
                        f"\t{args.sample_id}\n")
                    seen_header = True

                # handle data lines
                fields = line.strip().split('\t')
                if len(fields) == 1:
                    # this is a bad VCF file and the cols needs to be reassembled
                    if in_record:
                        # first col is the 'x' garbage, only append to
                        # current_record when in_record
                        current_record.append(line.strip())
                    else:
                        sys.stderr.write(
                            "[fix_vcf] Malformed QDNASeq VCF entry starting on line"
                            f"{line_i+1}, will squash lines and fix.\n")
                        in_record = True
                else:
                    # must have hit a good VCF line
                    if in_record:
                        # insert existing garbage before writing out good stuff
                        result.write('\t'.join(current_record) + '\n')
                        in_record = False
                        current_record = []

                    # pass through good lines
                    result.write(line)
        if in_record:
            result.write('\t'.join(current_record) + '\n')


if __name__ == "__main__":
    parser = argparser()
    args = parser.parse_args()
    main(args)
