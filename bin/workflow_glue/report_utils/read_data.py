"""Read data for report."""
import os

import pandas as pd
from pandas.api import types as pd_types

from .common import CATEGORICAL, CHROMOSOMES  # noqa: ABS101


def bamstats(stats_dir):
    """Read bamstats per-read stats.

    :param stats_dir: directory with bamstats per-read stats output files
    :return: `pd.DataFrame` with bamstats per-read stats data
    """
    # these files can be quite large; so only keep relevant columns and store the string
    # columns as categorical
    relevant_stats_cols_dtypes = {
        "name": str,
        "sample_name": CATEGORICAL,
        "ref": CATEGORICAL,
        "coverage": float,
        "ref_coverage": float,
        "read_length": int,
        "mean_quality": float,
        "acc": float,
    }
    input_files = os.listdir(stats_dir)
    dfs = []
    for fname in input_files:
        try:
            df = pd.read_csv(
                f"{stats_dir}/{fname}",
                sep="\t",
                index_col=0,
                usecols=relevant_stats_cols_dtypes,
                dtype=relevant_stats_cols_dtypes,
            )
        except pd.errors.EmptyDataError:
            df = pd.DataFrame(columns=relevant_stats_cols_dtypes)
        dfs.append(df)
    # check that we actually got any reads
    if not dfs:
        return pd.DataFrame(columns=relevant_stats_cols_dtypes)
    # pd.concat will revert the `dtype` of categorical columns back to `object` if the
    # concatenated columns don't contain the same set of items --> we use
    # `union_categoricals` to avoid that
    cat_cols = [k for k, v in relevant_stats_cols_dtypes.items() if v == CATEGORICAL]
    for col in cat_cols:
        uc = pd_types.union_categoricals([df[col] for df in dfs], sort_categories=True)
        for df in dfs:
            df[col] = pd.Categorical(df[col], categories=uc.categories, ordered=True)
    return pd.concat(dfs)


def flagstat(flagstat_dir):
    """Read bamstats per-file stats.

    :param stats_dir: directory with bamstats per-file stats output files
    :return: `pd.DataFrame` with bamstats per-file stats data
    """
    relevant_stats_cols_dtypes = {
        "ref": CATEGORICAL,
        "sample_name": CATEGORICAL,
        "total": int,
        "primary": int,
        "secondary": int,
        "supplementary": int,
        "unmapped": int,
        "qcfail": int,
        "duplicate": int,
    }
    input_files = os.listdir(flagstat_dir)
    dfs = []
    for fname in input_files:
        try:
            df = pd.read_csv(f"{flagstat_dir}/{fname}", sep="\t")
        except pd.errors.EmptyDataError:
            df = pd.DataFrame()
        dfs.append(df)
    if not dfs:
        return pd.DataFrame(columns=relevant_stats_cols_dtypes)
    return pd.concat(dfs)


def depths(depths_dir):
    """Read depth data.

    :param depths_dir: directory with `mosdepth` output files
    :param sample_names: list of sample names for which we found reads
    :return: `pd.DataFrame` with columns `["chrom", "start", "end", "depth",
        "sample_name"]`
    """
    dfs = []
    input_files = os.listdir(depths_dir)
    relevant_stats_cols_dtypes = {
        "chrom": CATEGORICAL,
        "start": int,
        "end": int,
        "depth": float,
    }
    for fname in input_files:
        sample_name = fname.split('.')[0]
        try:
            df = pd.read_csv(
                f"{depths_dir}/{fname}",
                sep="\t",
                names=relevant_stats_cols_dtypes,
                dtype=relevant_stats_cols_dtypes
            )
            df = df.eval(f'sample_name = "{sample_name}"')
            df = df.loc[df['chrom'].isin(CHROMOSOMES)]
        except pd.errors.EmptyDataError:
            df = pd.DataFrame(
                columns=relevant_stats_cols_dtypes.update({'sample_name': CATEGORICAL}))
        dfs.append(df)
    if not dfs:
        return pd.DataFrame(
            columns=relevant_stats_cols_dtypes.update({'sample_name': CATEGORICAL}))
    return pd.concat(dfs).astype({"sample_name": CATEGORICAL, "chrom": CATEGORICAL})
