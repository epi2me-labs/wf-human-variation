#!/usr/bin/env python
"""Plot CNVs."""

import base64

from dominate.tags import h6, img, p, span, table, tbody, td, th, thead, tr
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports.labs import LabsReport
from ezcharts.components.theme import LAB_head_resources
from ezcharts.layout.snippets import DataTable
from ezcharts.layout.snippets import Stats
from ezcharts.plots import Plot, util
from ezcharts.plots.distribution import histplot
import ezcharts.plots.ideogram as ideo
import pandas as pd
from pkg_resources import resource_filename
from .util import wf_parser  # noqa: ABS101

# Setup simple globals
WORKFLOW_NAME = 'wf-cnv'
REPORT_TITLE = f'{WORKFLOW_NAME} report'
Colors = util.Colors

GENETIC_SEX = {
    'XX': {
        'X': 0,
        'Y': -2},
    'XY': {
        'X': -1,
        'Y': -1},
    'XXX': {
        'X': 1,
        'Y': -2},
    'X0': {
        'X': -1,
        'Y': -2},
    'XXY': {
        'X': 0,
        'Y': -1}
}

KNOWN_CHROMOSOME_CONFIGURATIONS = {
    'Trisomy 21': {
        '1': 0,
        '2': 0,
        '3': 0,
        '4': 0,
        '5': 0,
        '6': 0,
        '7': 0,
        '8': 0,
        '9': 0,
        '10': 0,
        '11': 0,
        '12': 0,
        '13': 0,
        '14': 0,
        '15': 0,
        '16': 0,
        '17': 0,
        '18': 0,
        '19': 0,
        '20': 0,
        '21': 1,
        '22': 0,
    }
}


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("cnv_plot")
    parser.add_argument(
        '-q', '--qdna_seq', required=True, dest="qdnaseq_results",
        help="Output from running QDNAseq")
    # parser.add_argument('-s', '--sample', required=True, dest="sample",
    # help="Sample name")
    parser.add_argument(
        '-r', '--read_stats', required=True, dest="read_stats",
        help="Output from bamstats")
    parser.add_argument(
        '--bin_size', type=int,
        help="workflow paramters"
    )
    parser.add_argument(
        '--genome',
        help="genome version"
    )
    parser.add_argument(
        '--sample_id',
        help="sample_id"
    )
    parser.add_argument(
        '--noise_plot',
        help="QDNAseq noise plot"
    )
    parser.add_argument(
        '--isobar_plot',
        help="QDNAseq isobar plot"
    )
    parser.add_argument(
        '--params',
        help="workflow parameters"
    )
    parser.add_argument(
        '--versions',
        help="workflow versions"
    )
    parser.add_argument(
        '-o', '--output', required=True, dest="output_report",
        help="Output report")

    return parser


def process_qdnaseq_plots(plot):
    """Take qdnaseq PNG plot, return encoded string."""
    with open(plot, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read())

    return encoded_string


def process_qdna(qdnaseq_results):
    """Read qdnaseq results."""
    d = pd.read_csv(
        qdnaseq_results, sep="\t", skiprows=1,
        names=['chr', 'start', 'stop', 'bin', 'reads', 'strand', 'segs'])

    # get the most frequent seg count for each chromosome
    aggregate_segs = d.groupby(['chr'])['segs'].agg(pd.Series.mode)
    dict_of_calls = aggregate_segs.to_dict()

    return d, dict_of_calls


def process_fastcat(fastcat_file):
    """Load fastcat results into dataframe."""
    d = pd.read_csv(
        fastcat_file, sep="\t", header=0)

    return d


def process_chr_sizes(genome):
    """Load the chromosome sizes."""
    s = f"data/reference/{genome}/{genome}.chrom.sizes.gz"
    sizes = pd.read_csv(
        resource_filename(
            'ezcharts',
            s),
        sep="\t",
        header=None, names=['chr', 'length'])

    sizes['percent'] = sizes['length']/sizes['length'].sum()*88

    return sizes


def process_images(images):
    """Take a list of images and make into base64."""
    result = []
    for image in images:
        with open(image, "rb") as image_file:
            encoded_string = base64.b64encode(
                image_file.read()).decode('ascii')
            result.append(encoded_string)
    return result


def call_genetic_sex(chr_calls):
    """Call genetic sex from QDNAseq results."""
    genetic_sex = 'Undetermined'
    for call, items in GENETIC_SEX.items():
        if items['X'] == chr_calls['X']:
            if items['Y'] == chr_calls['Y']:
                genetic_sex = call

    return genetic_sex


def read_length_plot(read_lengths):
    """Make a histogram of read lengths."""
    histogram_data = read_lengths['read_length'].values

    plt = histplot(data=histogram_data, binwidth=100)

    xaxis = {
        'name': "Read length",
        'nameGap': '30',
        'nameTextStyle': {'fontSize': '16', 'fontStyle': 'bold'},
        'nameLocation': 'middle'
    }

    yaxis = {
        'name': "Number of reads",
        'nameGap': '60',
        'nameTextStyle': {'fontSize': '16', 'fontStyle': 'bold'},
        'nameLocation': 'middle'
    }

    plt.xAxis = xaxis
    plt.yAxis = yaxis

    return plt


def copy_number_plot(
        df_smoothed_read_counts,
        chr_sizes,
        style):
    """Create example plot."""
    # We need to set up a lot of lists!!
    plt = Plot()
    grid = list()
    xaxis = list()
    yaxis = list()
    left = 4
    left_increment = 0.2
    master_dataset_index = 0
    show = True
    # type = 'slider'
    map_cnv = list()
    # datazooms = list()
    xindexes = list()
    chromosomes = df_smoothed_read_counts.chr.unique()

    for count, chromosome in enumerate(chromosomes):
        # to do - grid sizes need to represent chr sizes
        width = chr_sizes.loc[
            chr_sizes['chr'] == f"chr{chromosome}"]["percent"]

        grid.append({
            'left': f'{left}%',
            'top': '3.5%',
            'width': f'{float(width)}%',
            'height': '85%',
            'show': True,
            'borderColor': Colors.grey90 if (count % 2) == 0 else Colors.white,
            'backgroundColor': Colors.grey90 if (
                count % 2) == 0 else Colors.white})

        left += float(width) + left_increment

        xaxis.append({
            'name': f"chr{chromosome}",
            'nameRotate': '45',
            'nameLocation': 'middle',
            'splitLine': {'show': False},
            'type': 'value',
            'gridIndex': count,
            'axisTick': {'show': False},
            'axisLabel': {'show': False},
            'splitNumber': 0,
            'show': True})

        xindexes.append(count)

        yaxis.append({
            'name': "Log2 Ratio",
            'nameLocation': 'middle',
            'gridIndex': count,
            'min': -2.1,
            'max': 2.1,
            'show': show,
            'splitLine': {'show': False},
            'type': 'value',
            'position': 'left',
            'axisLine': {'onZero': False}})
        show = False

        plt.add_series(dict(
            symbolSize=2,
            type='scatter',
            itemStyle={'color': Colors.grey70},
            xAxisIndex=count,
            yAxisIndex=count,
            datasetIndex=master_dataset_index
        ))

        chromsome_cnv_data = df_smoothed_read_counts.loc[
            df_smoothed_read_counts['chr'] == chromosome]

        data = []

        for index, row in chromsome_cnv_data.iterrows():
            data.append(
                [(row['start']+row['stop'])/2, row['reads']])

        # add track dataset
        plt.add_dataset(dict(source=data))

        # increment our dataset index
        master_dataset_index += 1

        #########
        # Plot the segments
        #########

        plt.add_series(dict(
            symbolSize=4,
            type='scatter',
            symbol='rect',
            itemStyle={'color': Colors.not_black},
            xAxisIndex=count,
            yAxisIndex=count,
            datasetIndex=master_dataset_index
        ))

        map_cnv.append(dict(
            left='right',
            show=False,
            min=-1,
            max=1,
            seriesIndex=master_dataset_index,
            inRange={
                'color': [Colors.cerulean, Colors.not_black, Colors.cinnabar]},
            text=['>0', '<0'],
            calculable=False))

        chromsome_seg_data = df_smoothed_read_counts.loc[
            df_smoothed_read_counts['chr'] == chromosome]

        data = []

        for index, row in chromsome_seg_data.iterrows():
            if row['reads'] != -1022:
                data.append(
                    [row['start'], row['segs']])

        plt.add_dataset(dict(source=data))

        # increment our dataset index
        master_dataset_index += 1

    plt.visualMap = map_cnv
    plt.grid = grid
    plt.xAxis = xaxis
    plt.yAxis = yaxis

    return plt


def make_report(
        read_lengths,
        chr_sizes,
        df_smoothed_read_counts,
        chr_calls,
        genetic_sex,
        bin_size,
        params,
        versions,
        genome,
        images,
        args):
    """Layout all the compeonets of the report."""
    report = LabsReport(
        f"{args.sample_id} | {REPORT_TITLE}", WORKFLOW_NAME, params, versions,
        head_resources=[*LAB_head_resources])

    with report.main_content:
        Stats(
            columns=3,
            items=[
                (genetic_sex, 'Chromosomal Sex'),
                (
                    f"{int(read_lengths['read_length'].median())}bp",
                    'Median Read Length'),
                (f'{"{:,}".format(len(read_lengths.index))}', 'Total Reads')
            ])

    with report.add_section(
            'Chromosome Copy Summary', 'Summary'):

        conversion = {0: -2, 1: -1, 2: 0, 3: 1, 4: 2}
        with table(cls="table"):
            with thead():
                for copies, cnv_call in conversion.items():
                    th(f"{copies} copies")
            with tbody():
                with tr():
                    for copies, cnv_call in conversion.items():
                        cell = td()
                        count = 0
                        for chrom, call in chr_calls.items():
                            if call == cnv_call:
                                count += 1
                                if call > 0:
                                    cell.add(
                                        span(
                                            chrom,
                                            cls="badge bg-danger"))
                                elif call < 0:
                                    cell.add(
                                        span(
                                            chrom,
                                            cls="badge bg-primary"))
                                else:
                                    cell.add(
                                        span(
                                            chrom,
                                            cls="badge bg-light text-dark"))
                        if count == 0:
                            cell.add("-")

    with report.add_section(
            'Ideoplot', 'Ideoplot'):

        cnv_data = pd.read_csv(
            args.qdnaseq_results,
            sep="\t",
            skiprows=1,
            header=None,
            names=['chr', 'start', 'end', 'name', 'value', 'strand', 'call'])

        cnv_colors = {
            -2: util.Colors.cerulean,
            -1: util.Colors.cerulean,
            0: util.Colors.grey40,
            1: util.Colors.cinnabar,
            2: util.Colors.cinnabar
        }

        cnv_data['color'] = cnv_data['call'].map(cnv_colors)
        cnv_data.loc[cnv_data.value == -1022, 'color'] = util.Colors.sandstorm

        cnv_data = cnv_data.loc[cnv_data['call'] != 0]

        track_data = pd.read_csv(
            args.qdnaseq_results,
            sep="\t",
            skiprows=1,
            header=None,
            names=['chr', 'start', 'end', 'name', 'value', 'strand', 'call'])

        track_data['color'] = util.Colors.grey80

        band_types = ['acen', 'gvar', 'stalk']
        bands_df = pd.read_csv(
            resource_filename(
                'ezcharts',
                f"data/reference/{genome}/cytoBand.txt.gz"),
            sep="\t",
            header=None)

        bands_df = bands_df.set_axis([
            'chr',
            'start',
            'end',
            'band',
            'type1'
        ], axis=1)

        bands_df = bands_df[bands_df['type1'].isin(band_types)]
        bands_df['strand'] = 'NA'
        bands_df['type'] = bands_df['type1']
        bands_df['color'] = util.Colors.not_black
        bands_df['chr'] = bands_df['chr'].str.replace('chr', '')

        EZChart(
            ideo.ideogram(
                blocks=[bands_df, cnv_data],
                track=track_data,
                genome=genome),
            height='1000px', width='100%', theme='epi2melabs')
        p("""Red: Increased copy number detected""")
        p("""Blue: Decreased copy number detected""")
        p("""Yellow: QDNASeq no call""")
        p("""Black: 'acen', 'gvar', or 'stalk'.
            Read counts of < -2 converted to 2""")

    with report.add_section(
            'Detailed CNV Results', 'Details'):

        p("QDNASeq Smoothed Read Counts")
        p(f"Bin size {bin_size}Kbp, No call not displayed")
        EZChart(
            copy_number_plot(
                df_smoothed_read_counts, chr_sizes, "scatter"),
            'epi2melabs')

        data_table = DataTable(
            headers=[
                'chr',
                'start',
                'end',
                'genes',
                'read count',
                'call'])

        calls = df_smoothed_read_counts.loc[
            df_smoothed_read_counts['segs'] != 0]

        for index, row in calls.iterrows():
            ensembl = "https://rest.ensembl.org/overlap/region/human/"
            region = f"{row['chr']}:{row['start']}-{row['stop']}"
            options = "?feature=gene;content-type=text/x-gff3"
            genes_in_region = f"""
                <a href=\"{ensembl}{region}{options}\">Genes</a>"""
            data_table.add_row(
                title=None,
                columns=[
                    f"chr{row['chr']}",
                    row['start'],
                    row['stop'],
                    genes_in_region,
                    row['reads'],
                    row['segs']])

    with report.add_section('Quality Control Data', 'QC'):
        h6("Histogram of Read Lengths")
        EZChart(read_length_plot(read_lengths), 'epi2melabs')
        p("""N.B. The read length affects the results of QDNASeq analysis.
            In future versions we will provide preset parameters based on
            detected read length""")

        h6("QDNASeq Copy Number QC Plots")
        with table():
            with tbody():
                with tr():
                    for image in images:
                        td(img(src=f"data:image/png;base64,{image}"))
                with tr():
                    td("""Median read counts per genomic bin shown as a
                        function of GC content and mappability""")
                    td("""Read count relationship between sequence depth
                        and noise""")

    return report


def main(args):
    """Run the entry point."""
    df_smoothed_read_counts, chr_calls = process_qdna(args.qdnaseq_results)
    genetic_sex = call_genetic_sex(chr_calls)
    read_lengths = process_fastcat(args.read_stats)
    chr_sizes = process_chr_sizes(args.genome)
    images = process_images([args.isobar_plot, args.noise_plot])

    report = make_report(
        read_lengths,
        chr_sizes,
        df_smoothed_read_counts,
        chr_calls,
        genetic_sex,
        args.bin_size,
        args.params,
        args.versions,
        args.genome,
        images,
        args)

    report.write(args.output_report)
