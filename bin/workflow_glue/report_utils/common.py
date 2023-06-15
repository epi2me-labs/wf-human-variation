"""Common variables."""

from ezcharts.plots import util
from ezcharts.plots.distribution import histplot
import numpy as np
from pandas.api import types as pd_types

CATEGORICAL = pd_types.CategoricalDtype(ordered=True)
CHROMOSOMES = {str(i): i for i in range(1, 23)}
CHROMOSOMES.update({f"chr{i}": int(i) for i in CHROMOSOMES})
CHROMOSOMES.update({'X': 23, 'Y': 24, 'chrX': 23, 'chrY': 24})
THEME = "epi2melabs"
CLINVAR_BASE = "https://www.ncbi.nlm.nih.gov/clinvar/variation/"
NCBI_GENE_BASE = "https://www.ncbi.nlm.nih.gov/gene/"


def hist_plot(
        df, col, title, xaxis='', yaxis='', rounding=0,
        color=None, binwidth=None, max_y=None, max_x=None):
    """Make a histogram of given parameter."""
    histogram_data = df[col].values
    plt = histplot(data=histogram_data, binwidth=binwidth, color=color)
    meanv = np.mean(histogram_data).round(rounding)
    medianv = np.median(histogram_data).round(rounding)

    # Define the title
    plt.title = dict(
        text=title,
        subtext=(
            f"Mean: {meanv}. "
            f"Median: {medianv}. "
        ),
    )

    # Add median value (Thanks Julian!)
    plt.add_series(
        dict(
            type="line",
            name="Mean",
            data=[dict(value=[meanv, 0]), dict(value=[meanv, max_y])],
            itemStyle=(dict(color=util.Colors.sandstorm))
        )
    )
    # Add median value (Thanks Julian!)
    plt.add_series(
        dict(
            type="line",
            name="Median",
            data=[dict(value=[medianv, 0]), dict(value=[medianv, max_y])],
            itemStyle=(dict(color=util.Colors.fandango))
        )
    )
    # Customize X-axis
    xaxis = {
        'name': xaxis,
        'nameGap': '30',
        'nameTextStyle': {'fontSize': '16', 'fontStyle': 'bold'},
        'nameLocation': 'middle'
    }
    # Customize Y-axis
    yaxis = {
        'name': yaxis,
        'nameGap': '15',
        'nameTextStyle': {'fontSize': '16', 'fontStyle': 'bold'},
        'nameLocation': 'middle'
    }
    # Update max Y value if max_y is passed
    if max_y:
        yaxis.update(max=max_y)
    # Update max X value if max_x is passed
    if max_x:
        xaxis.update(max=max_x)
    # Update axis parameters
    plt.xAxis = xaxis
    plt.yAxis = yaxis

    # Print series
    for s in plt.series:
        if s.type == 'line':
            s.symbolSize = 0
    return plt
