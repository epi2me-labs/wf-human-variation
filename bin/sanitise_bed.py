#!/usr/bin/env python
"""Compare chromosome labels in BED and REF, and error if mismatched."""

import argparse
import os
import sys

from natsort import natsort_keygen
import pandas as pd
from pyfaidx import Fasta


def read_bed_file(bed_file):
    """Read BED file, return df."""
    bed_data = pd.read_csv(bed_file, delim_whitespace=True, header=None)
    return bed_data


def get_chrs_from_fasta(fasta_file):
    """Read reference FASTA file, return the chromosome labels."""
    fasta = Fasta(fasta_file)
    return list(fasta.keys())


def write_updated_bed_file(updated_bed_data, output_bed_file):
    """Write updated BED data to a new BED file."""
    updated_bed_data.to_csv(output_bed_file, sep='\t', index=False, header=False)


def main():
    """Run entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ref",
        help="Reference FASTA file",
        required=True
    ),
    parser.add_argument(
        "--bed",
        help="BED file",
        required=True
    )
    parser.add_argument(
        "--bed_out",
        help="Sanitised BED file",
        required=True
    )
    args = parser.parse_args()

    bed_data = read_bed_file(args.bed)
    fasta_chrom_labels = get_chrs_from_fasta(args.ref)

    # unique list of BED chrs, converted to str so they can be compared later
    unique_bed_chromosomes = bed_data[0].unique().astype(str).tolist()

    # check all bed chromosomes are in fasta chromsome list
    if set(unique_bed_chromosomes).issubset(set(fasta_chrom_labels)):
        # sort and write to new file with correct whitespace:
        bed_data.sort_values(
            by=[bed_data.columns[0], bed_data.columns[1]],
            key=natsort_keygen(),
            inplace=True
        )
        write_updated_bed_file(bed_data, args.bed_out)
    # discrepancy between ref and bed chromosomes
    else:
        sys.stderr.write(
            "There is a discrepancy in chromosome labelling between the "
            "reference FASTA and BED file. Perhaps you have forgotten to "
            "add or remove the chr prefix?\n")
        sys.exit(os.EX_DATAERR)


if __name__ == '__main__':
    main()
