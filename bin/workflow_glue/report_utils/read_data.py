"""Read data for report."""
import os

import pandas as pd
from pandas.api import types as pd_types
from pysam import VariantFile

from .common import CATEGORICAL, CHROMOSOMES  # noqa: ABS101


def fasta_idx(faidx):
    """Read faidx for the reference fasta."""
    # these files can be quite large; so only keep relevant columns and store the string
    # columns as categorical
    relevant_stats_cols_dtypes = {
        "chrom": CATEGORICAL,
        "length": int,
        "offset1": int,
        "offset2": int,
        "offset3": int,
    }
    try:
        df = pd.read_csv(
            faidx,
            sep="\t",
            names=relevant_stats_cols_dtypes,
            dtype=relevant_stats_cols_dtypes,
            usecols=['chrom', 'length']
        )
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=relevant_stats_cols_dtypes)
    # check that we actually got any reads
    return df


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


def make_breaks(minval, maxval, winsize):
    """Create intervals inclusive of last value."""
    # Create the intervals
    breaks = list(range(minval, maxval, winsize))
    # Check that the max value is the chr length
    if breaks[-1] > maxval:
        breaks[-1] = maxval
    if breaks[-1] < maxval:
        breaks += [maxval]
    return breaks


def add_missing_windows(intervals, faidx, winsize=25000):
    """Add missing windows to depth dataframe."""
    relevant_stats_cols_dtypes = {
        "chrom": CATEGORICAL,
        "start": int,
        "end": int,
        "depth": float,
    }
    # Get unique chromosomes
    chrs = intervals.chrom.unique().tolist()
    # New intervals
    final_intervals = []
    # Add missing windows for each chromosome
    for chr_id in chrs:
        # Get chromosome entries
        chr_entry = intervals.loc[intervals['chrom'] == chr_id].reset_index(drop=True)
        # First window start
        minval = chr_entry.start.min()
        # Final window end
        maxval = chr_entry.end.max()
        # Chromosome length
        chr_len = faidx[faidx['chrom'] == chr_id].length.max()
        # If the minimum value is !=0, add intermediate windows
        if minval != 0:
            # Create the breaks for the intervals
            breaks = make_breaks(0, minval, winsize)
            # Create vectors to populate the DF
            chr_vec = [chr_id for i in range(0, len(breaks) - 1)]
            depth_vals = [0 for i in range(0, len(breaks) - 1)]
            # Add heading    intervals
            final_intervals.append(
                pd.DataFrame(data={
                    'chrom': chr_vec,
                    'start': breaks[0:-1],
                    'end': breaks[1:],
                    'depth': depth_vals}))
        # Append precomputed intervals, checking for breaks
        # in between regions.
        for idx, region in chr_entry.iterrows():
            # To DF
            region = region.to_frame().T
            # Add the first window as default
            if idx == 0:
                final_intervals.append(region)
                continue
            # If the new window start is not the end of the previous,
            # then create new regions
            if region.start.min() != final_intervals[-1].end.max():
                # Create intervals
                breaks = make_breaks(
                    final_intervals[-1].end.max(),
                    region.start.min(),
                    winsize)
                # Create vectors to populate the DF
                chr_vec = [chr_id] * (len(breaks) - 1)
                depth_vals = [0] * (len(breaks) - 1)
                # Add heading    intervals
                final_intervals.append(
                    pd.DataFrame(data={
                        'chrom': chr_vec,
                        'start': breaks[0:-1],
                        'end': breaks[1:],
                        'depth': depth_vals}))
                final_intervals.append(region)
            else:
                final_intervals.append(region)
        # If the max value is less than the chromosome length, add these too
        if maxval < chr_len:
            # Create the breaks for the intervals
            breaks = make_breaks(maxval, chr_len, winsize)
            # Create vectors to populate the DF
            chr_vec = [chr_id] * (len(breaks) - 1)
            depths = [0] * (len(breaks) - 1)
            # Add tailing intervals
            final_intervals.append(
                pd.DataFrame(data={
                    'chrom': chr_vec,
                    'start': breaks[0:-1],
                    'end': breaks[1:],
                    'depth': depths}))
    return pd.concat(final_intervals).astype(
        relevant_stats_cols_dtypes).reset_index(drop=True)


def depths(depths_dir, faidx, winsize):
    """Read depth data.

    :param depths_dir: directory with `mosdepth` output files
    :param faidx: dataframe from the fasta fai index
    :winsize: size of the missing windows to add
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
            df = add_missing_windows(df, faidx, winsize)
            df = df.eval(f'sample_name = "{sample_name}"')
            df = df.loc[df['chrom'].isin(CHROMOSOMES)]
        except pd.errors.EmptyDataError:
            df = pd.DataFrame(
                columns=relevant_stats_cols_dtypes.update({
                    'sample_name': CATEGORICAL,
                    'length': int
                    }))
        dfs.append(df)
    if not dfs:
        return pd.DataFrame(
            columns=relevant_stats_cols_dtypes.update({
                'sample_name': CATEGORICAL,
                'length': int}))
    return pd.concat(dfs).astype({"sample_name": CATEGORICAL, "chrom": CATEGORICAL})


def parse_vcf(fname):
    """Parse input VCF using pysam."""
    vcf_file = VariantFile(fname)
    all_variants = vcf_file.fetch()
    return (all_variants)


def parse_vcf_for_size(fname):
    """Read VCF into a dataframe."""
    try:
        df = pd.read_csv(fname, comment='#', nrows=10, header=None, delimiter='\t')
        return False, df
    except pd.errors.EmptyDataError:
        df = pd.DataFrame()
        return True, df
