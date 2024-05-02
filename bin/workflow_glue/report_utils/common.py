"""Common variables."""

import os

from ezcharts.plots import util
from ezcharts.plots.distribution import histplot
import numpy as np
import pandas as pd
from pandas.api import types as pd_types

CATEGORICAL = pd_types.CategoricalDtype(ordered=True)
CHROMOSOMES = {str(i): i for i in range(1, 23)}
CHROMOSOMES.update({f"chr{i}": int(i) for i in CHROMOSOMES})
CHROMOSOMES.update({'X': 23, 'Y': 24, 'chrX': 23, 'chrY': 24})
THEME = "epi2melabs"


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


def compute_n50(data, x='start', y='count'):
    """Automatic detection of the data type for N50."""
    if isinstance(data, pd.DataFrame):
        n50_value = n50_hist(data, x=x, y=y)
    elif isinstance(data, np.ndarray):
        n50_value = n50_array(data)
    elif isinstance(data, list) or isinstance(data, tuple):
        n50_value = n50_array(np.array(data))
    else:
        raise ValueError(f"Unsupported data type {type(data)}")
    return n50_value


def n50_array(lengths):
    """Compute read N50.

    :param length: numpy vector of lengths
    """
    # Sort the read lengths
    sorted_l = np.sort(lengths)[::-1]
    # Generate cumsum
    cumsum = np.cumsum(sorted_l)
    # Get lowest cumulative value >= (total_length/2)
    n50 = sorted_l[np.searchsorted(cumsum, cumsum[-1]/2)]
    return n50


def n50_hist(length_hist, x='start', y='count'):
    """Compute read N50 from histogram data."""
    # Create the vector from the histogram data
    cumsum = np.cumsum(length_hist[y].values * length_hist[x].values)
    # Get lowest cumulative value >= (total_length/2)
    n50 = length_hist.iloc[np.searchsorted(cumsum, cumsum[-1]/2)].start
    return n50


def sum_hists(mapped_hist, unmapped_hist):
    """Sum two histogram dataframes based on the intervals."""
    return pd.concat(
        [mapped_hist, unmapped_hist]
    )\
        .reset_index(drop=True)\
        .sort_values(by='start')\
        .groupby(['start', 'end'])\
        .sum().reset_index()


def load_hists(hists_dir, dtype):
    """Load and combine histogram data."""
    # Data type
    dt = float
    if "length" in dtype:
        dt = int

    # Load mapped
    hist_map = pd.read_csv(
        os.path.join(hists_dir, f"{dtype}.hist"),
        sep="\t",
        names=["start", "end", "count"],
        dtype={"start": dt, "end": dt, "count": int},
    )
    # Try loading unmapped
    if dtype in ['length', 'quality']:
        hist_umap = pd.read_csv(
            os.path.join(hists_dir, f"{dtype}.unmap.hist"),
            sep="\t",
            names=["start", "end", "count"],
            dtype={"start": dt, "end": dt, "count": int},
        )
    else:
        hist_umap = pd.DataFrame(
            names=["start", "end", "count"],
            dtype={"start": dt, "end": dt, "count": int},
        )
    # Sum mapped and unmapped histograms.
    hist = sum_hists(hist_map, hist_umap)
    return hist, hist_map, hist_umap
