#!/usr/bin/env python
"""Fix malformed QDNAseq VCFs by reformatting files with a single CNV call."""

from .util import wf_parser  # noqa: ABS101


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("get_genome")
    parser.add_argument(
        '--vcf', required=True, dest="vcf",
        help="Malformed VCF")
    parser.add_argument(
        '--fixed_vcf', required=True, dest="fixed_vcf",
        help="Fixed VCF")
    parser.add_argument(
        '--sample_id', required=True, dest="sample_id")
    return parser


def main(args):
    """Run entry point."""
    vcf_record = []
    result = open(args.fixed_vcf, 'w')
    with open(args.vcf, 'r') as vcf_file:
        for line in vcf_file:
            line = line.rstrip()
            if line.startswith("#"):
                result.write(line+"\n")
            else:
                vcf_record.append(line)

    result.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
    result.write(args.sample_id+"\n")

    # Remove the rogue 'x' that appears in the .seg file when only one CNV is called
    vcf_record.pop(0)

    result.write('\t'.join(vcf_record))
