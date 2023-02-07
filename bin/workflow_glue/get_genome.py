#!/usr/bin/env python
"""Check BAM header and return genome build based on chromosome sizes."""

from .util import wf_parser  # noqa: ABS101


CHROMOSOME_SIZES = {
    'hg19': {
        'chr1': '249250621',
        'chr2': '243199373',
        'chr3': '198022430',
        'chr4': '191154276',
        'chr5': '180915260',
        'chr6': '171115067',
        'chr7': '159138663',
        'chr8': '146364022',
        'chr9': '141213431',
        'chr10': '135534747',
        'chr11': '135006516',
        'chr12': '133851895',
        'chr13': '115169878',
        'chr14': '107349540',
        'chr15': '102531392',
        'chr16': '90354753',
        'chr17': '81195210',
        'chr18': '78077248',
        'chr19': '59128983',
        'chr20': '63025520',
        'chr21': '48129895',
        'chr22': '51304566',
        'chrX': '155270560',
        'chrY': '59373566'},
    'hg38': {
        'chr1': '248956422',
        'chr2': '242193529',
        'chr3': '198295559',
        'chr4': '190214555',
        'chr5': '181538259',
        'chr6': '170805979',
        'chr7': '159345973',
        'chr8': '145138636',
        'chr9': '138394717',
        'chr10': '133797422',
        'chr11': '135086622',
        'chr12': '133275309',
        'chr13': '114364328',
        'chr14': '107043718',
        'chr15': '101991189',
        'chr16': '90338345',
        'chr17': '83257441',
        'chr18': '80373285',
        'chr19': '58617616',
        'chr20': '64444167',
        'chr21': '46709983',
        'chr22': '50818468',
        'chrX': '156040895',
        'chrY': '57227415'}
}


ALLOWED_CHR = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
    "chr10", "chr11", "chr12", "chr11", "chr12", "chr13", "chr14", "chr15",
    "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX",
    "chrY"
]


def chromosome_sizes(extracted_sizes):
    """Get dictionary of chromosomes and sizes from samtools index."""
    fasta_sizes = {}
    with open(extracted_sizes, 'r') as csvfile:
        for line in csvfile:
            line = line.rstrip()
            cols = line.split('\t')
            if cols[0] in ALLOWED_CHR:
                fasta_sizes[cols[0]] = cols[1]

    return fasta_sizes


def get_genome(sizes):
    """Get genome based on chromosome sizes."""
    genome = ""
    for i in CHROMOSOME_SIZES.keys():
        if sizes.items() <= CHROMOSOME_SIZES[i].items():
            genome = i
    return genome


def main(args):
    """Run entry point."""
    all_sizes = chromosome_sizes(args.chr_counts)

    genome = get_genome(all_sizes)
    result = open(args.output, 'w')
    result.write(genome)


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("get_genome")
    parser.add_argument(
        '--chr_counts', required=True, dest="chr_counts",
        help="Output from samtools faidx")
    parser.add_argument(
        '-o', '--output', dest="output",
        help="Output genome")
    return parser
