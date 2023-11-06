#!/usr/bin/env python
"""Generate data input files for STR content plot."""

import json
import re

import pandas as pd
import pysam

from .util import wf_parser  # noqa: ABS101


def argparser():
    """Return an arg parser object from arguments."""
    parser = wf_parser("generate_str_content")

    parser.add_argument(
        '--straglr',
        help="STRAGLR tsv",
        required=True
    )
    parser.add_argument(
        '--stranger',
        help="STRANGER tsv",
        required=True
    )
    parser.add_argument(
        '--chr',
        help="BAM Chr",
        required=True
    )
    parser.add_argument(
        '--repeat_bed',
        help="BED file of STR repeats",
        required=True
    )
    parser.add_argument(
        '--str_reads_bam',
        help="BAM file with STR reads",
        required=True
    )

    return parser


def truncate_interruption(interruption_seq):
    """Shorten long interruption sequences for display on Bokeh HoverTool."""
    if len(interruption_seq) > 20:
        return interruption_seq[:21] + "..."
    else:
        return interruption_seq


def bed_ru_merge(merged_tsv, repeats_bed):
    """Merge RU's from STR repeats BED into straglr/stranger merged TSV."""
    bed_columns = [
        'bed_chr', 'bed_start', 'bed_end', 'bed_ru', 'bed_gene', 'bed_varid'
    ]
    remove = ['bed_start', 'bed_end', 'bed_gene']
    pd_repeats_bed = pd.read_csv(
        repeats_bed,
        sep="\t",
        header=None,
        names=bed_columns
    )
    pd_repeats_bed = pd_repeats_bed.drop(columns=remove)
    merged_tsv = pd.merge(
        merged_tsv,
        pd_repeats_bed,
        left_on='VARID',
        right_on='bed_varid'
    )
    merged_tsv = merged_tsv.drop(columns='bed_varid')
    return merged_tsv


def extract_sequences(bam, merged_tsv):
    """Extract STR sequences and related info from BAM and create JSON."""
    str_seq_dict = {}

    # Need to iterate through VARID's as
    # Straglr TSV contains duplicate read ID's supporting various VARID's
    unique_varids = merged_tsv['VARID'].unique()

    for varid in unique_varids:
        # Gathering supporting read_ids for the VARID
        relevant_read_ids = merged_tsv.loc[
            merged_tsv['VARID'] == varid, 'read'
        ].values

        input_bam = pysam.AlignmentFile(bam, "rb")

        for record in input_bam:
            if record.query_name in relevant_read_ids:
                read = merged_tsv['read'] == record.query_name
                merged_varid = merged_tsv[read & (merged_tsv['VARID'] == varid)]

                strand = merged_varid['strand'].values[0]
                chrom = merged_varid['bed_chr'].values[0]
                repeat_start = merged_varid['read_start'].values[0]
                size = merged_varid['size'].values[0]
                str_normal_max = int(merged_varid['STR_NORMAL_MAX'].values[0])
                str_pathologic_min = int(merged_varid['STR_PATHOLOGIC_MIN'].values[0])

                # RU from Repeats BED
                repeat_unit = merged_varid['bed_ru'].values[0]

                # Create STR identifier for plot title
                disease = merged_varid['Disease'].values[0]
                str_identifier = f"{disease} ({varid})"

                # Collate STR Summary info into dict
                if str_identifier not in str_seq_dict:
                    str_seq_dict[str_identifier] = {
                        'chr': chrom,
                        'repeat_unit': repeat_unit,
                        'VARID': varid,
                        'str_normal_max': str_normal_max,
                        'str_pathologic_min': str_pathologic_min,
                        'observed_reads': {}
                    }

                # Extract Haplotype
                if 'HP' in dict(record.tags):
                    haplotype = dict(record.tags)['HP']
                else:
                    haplotype = "None"

                # Extracting read and STR sequences
                # -1 as size includes repeat_start pos.
                repeat_end = (repeat_start + size) - 1
                read_sequence = record.query_sequence
                if read_sequence is not None:
                    if strand == "+":
                        str_sequence = read_sequence[
                            (repeat_start):(repeat_end + 1)
                        ]
                    if strand == "-":
                        reversed_repeat_start = len(read_sequence) - repeat_end
                        reversed_repeat_end = len(read_sequence) - repeat_start
                        str_sequence = read_sequence[
                            (reversed_repeat_start):(reversed_repeat_end + 1)
                        ]
                else:
                    continue

                # Detecting the RU's and interruptions
                repeat_unit_indexes = []
                interruption_indexes = []

                ru_regex = re.compile(fr'{repeat_unit}')
                for repeat_unit in ru_regex.finditer(str_sequence):
                    repeat_unit_indexes.append((
                        repeat_unit.start(),
                        repeat_unit.end(),
                        repeat_unit.group(),
                        repeat_unit.end() - repeat_unit.start()
                    ))

                # Finding interruptions at start of seq and between RU's
                previous_end = 0
                for start, end, seq, length in repeat_unit_indexes:
                    if start > previous_end:
                        interruption_indexes.append((
                            previous_end,  # start of interruption
                            start,  # start of RU/end of interruption
                            str_sequence[previous_end:start],  # interruption seq
                            start - previous_end  # len of interruption
                        ))
                    previous_end = end

                # Finding interruptions at end of sequence
                # If end of last RU is smaller than the seq len
                str_sequence_l = len(str_sequence)
                if previous_end < str_sequence_l:
                    interruption_indexes.append((
                        previous_end,
                        str_sequence_l,
                        str_sequence[previous_end:len(str_sequence)],
                        str_sequence_l - previous_end
                    ))

                str_seq_dict[str_identifier]['observed_reads'].update({
                    record.query_name: {
                        "str_sequence": str_sequence,
                        "str_seq_length": str_sequence_l,
                        "haplotype": haplotype,
                        "repeat_unit_indexes": repeat_unit_indexes,
                        "interruption_indexes": interruption_indexes
                    }
                })

        input_bam.close()

    str_seq_json = json.dumps(str_seq_dict, indent=4)
    return str_seq_json


def create_plot_input_files(str_seq_json):
    """Extract info relevant for plots from JSON and save as CSV."""
    data = json.loads(str_seq_json)

    for str_identifier, str_data in data.items():
        rows = []

        for read_id, read_details in str_data["observed_reads"].items():
            for seq_type, indexes in {
                "Repeat": read_details["repeat_unit_indexes"],
                "Interruption": read_details["interruption_indexes"]
            }.items():

                for index in indexes:
                    rows.append([
                        str_identifier,
                        str_data['chr'],
                        str_data['VARID'],
                        str_data['repeat_unit'],
                        str_data['str_normal_max'],
                        str_data['str_pathologic_min'],
                        read_id,
                        read_details["haplotype"],
                        read_details["str_seq_length"],
                        seq_type,
                        index[2],  # sequence
                        truncate_interruption(index[2]),
                        index[0],  # start pos
                        index[1],  # end pos
                        index[3]])  # length

        df = pd.DataFrame(rows, columns=[
            "str_identifier",
            "chromosome",
            "varid",
            "repeat_unit",
            "str_normal_max",
            "str_pathologic_min",
            "read_id",
            "haplotype",
            "str_seq_length",
            "type",
            "sequence",  # RU/interruption seq
            "truncated_seq",  # seq shortened for hover tool display
            "start",  # start pos in STR seq
            "end",  # end pos in STR seq
            "length"  # len of each RU/interruption
        ])

        str_identifier = str_identifier.replace(" (", "_")
        str_identifier = str_identifier.replace(")", "")
        str_identifier = str_identifier.replace(" ", "")
        df.to_csv(f"{str_identifier}_str-content.csv", index=False)


def main(args):
    """Run the entry point."""
    # Merging Straglr and Stranger Plot TSV's
    pd_straglr = pd.read_csv(args.straglr, sep="\t", header=1)
    pd_stranger = pd.read_csv(args.stranger, sep="\t", header=0)
    merged = pd.merge(pd_straglr, pd_stranger, left_on='start', right_on='POS')
    # Merging BED RU info
    merged = bed_ru_merge(merged, args.repeat_bed)
    # Create JSON with STR sequence info
    str_seq_json = extract_sequences(args.str_reads_bam, merged)
    # Create CSV's to be used as inputs for repeat content plots
    create_plot_input_files(str_seq_json)
