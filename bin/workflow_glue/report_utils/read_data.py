"""Read data for report."""
import os

from dominate.tags import a, span
from ezcharts.layout.util import isolate_context
from natsort import natsort_keygen
import pandas as pd
from pandas.api import types as pd_types
from pysam import VariantFile

from .common import CATEGORICAL, CHROMOSOMES  # noqa: ABS101
from .common import CLINVAR_BASE, NCBI_BASE  # noqa: ABS101


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


# vcf to pysam
def vcf_to_pysam(fname):
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


def clinvar_to_df(vcf_data):
    """Convert ClinVar VCF to sorted df."""
    vcf_file = VariantFile(vcf_data)
    all_variants = vcf_file.fetch()

    data = []
    columns = [
        "Chrom", "Pos", "Gene(s)", "ClinVar",
        "Significance", "Type", "Consequence", "HGVSc", "HGVSp"]
    for variant in all_variants:
        # clinvar significance - ignore if benign or likely benign
        significance = ", ".join(variant.info['CLNSIG'])
        if 'benign' in significance.lower():
            continue
        else:
            # NCBI gene URL
            try:
                all_ncbi_urls = []
                clinvar_gene_string = variant.info['GENEINFO']
                all_genes = clinvar_gene_string.split('|')
                for gene in all_genes:
                    gene_symbol, gene_id = gene.split(':')
                    with isolate_context():
                        ncbi_url = str(a(
                            gene_symbol,
                            href=f'{NCBI_BASE}{gene_id}'))
                        all_ncbi_urls.append(ncbi_url)
                ncbi_gene_url = ", ".join(all_ncbi_urls)
            except KeyError:
                ncbi_gene_url = "No affected genes found"

            # multiple ClinVar IDs possible, separated by ';'
            clinvar_url_list = []
            clinvar_id = variant.id
            all_clinvar_ids = clinvar_id.split(';')
            for each_id in all_clinvar_ids:
                with isolate_context():
                    clinvar_url = str(a(
                        each_id, href=f'{CLINVAR_BASE}{each_id}'))
                    clinvar_url_list.append(clinvar_url)
            clinvar_urls_for_report = ", ".join(clinvar_url_list)

            # tidy up significances
            significance = significance.replace("_", " ").replace(
                "|", ", ").capitalize()

            # tidy up variant types
            variant_type = variant.info['CLNVC']
            variant_type = variant_type.replace("_", " ").capitalize()
            if variant_type == 'Single nucleotide variant':
                variant_type = 'SNV'

            # tidy up consequences
            consequences = []
            try:
                all_consequences = variant.info['MC']
                for each_conseq in all_consequences:
                    ontology, consequence = each_conseq.split('|')
                    consequence = consequence.replace("_", " ").capitalize()
                    consequence = consequence.replace(
                        " prime utr", "' UTR")
                    consequences.append(consequence)
                consequences = ", ".join(consequences)
            except KeyError:
                consequences = "No consequences found"

            # HGVS p. and c. come from the SnpEff ANN field, we will take the first
            # available 'NM_' number
            hgvsc = ""
            hgvsp = ""
            all_variant_info = variant.info['ANN']
            for each_variant in all_variant_info:
                fields = each_variant.split('|')
                # find the first NM (because there are some XM records present)
                if fields[6].startswith("NM"):
                    hgvsc = f"{fields[6]}:{fields[9]}"
                    hgvsp = fields[10] if fields[10] != "" else "-"
                    break
                # if there are only XMs then display the first one of those instead
                else:
                    hgvsc = f"{fields[6]}:{fields[9]}"
                    hgvsp = fields[10] if fields[10] != "" else "-"

        record = (
            variant.chrom, variant.pos, ncbi_gene_url, clinvar_urls_for_report,
            significance, variant_type, consequences, hgvsc, hgvsp
        )
        data.append(record)

    df = pd.DataFrame(data, columns=columns)

    # get unique list of all significances in the df
    sample_significances = df['Significance'].unique()

    # desired order of ClinVar significances
    significance_order = [
        "Pathogenic", "Pathogenic, low penetrance",
        "Likely pathogenic", "Likely pathogenic, low penetrance",
        "Uncertain significance", "Conflicting interpretations of pathogenicity",
        "Risk factor", "Established risk allele", "Likely risk allele",
        "Uncertain risk allele", "Association", "Confers sensitivity", "Affects",
        "Protective", "Association not found", "Drug response", "Other",
        "Not provided"]

    """re-order variant significances found in the sample to the desired order
    also ensures there are no duplicate significances as this will cause a
    problem when setting 'Significance' to a categorical column later
    """
    reordered = [
        y for x in significance_order for y in sample_significances if y.startswith(x)]
    reordered_unique = sorted(set(reordered), key=reordered.index)

    # define custom sort for chromosomes
    def nat_sort_chromosome(col):
        if col.name == "Chrom":
            return col.apply(natsort_keygen())
        return col

    # sort the sample ClinVar results using the re-ordered significances
    df["Significance"] = pd.Categorical(
        df["Significance"], categories=reordered_unique)
    df.sort_values(
        by=["Significance", "Chrom", "Pos"], inplace=True, key=nat_sort_chromosome)

    return df


def format_clinvar_table(df):
    """Format ClinVar dataframe to add badges where required."""
    df_for_table = df.copy()
    # uncategorise the 'Significance' column so we can add badges
    df_for_table['Significance'] = df_for_table['Significance'].astype(
        df_for_table['Significance'].cat.categories.to_numpy().dtype)

    significance_badges = {
        "Pathogenic": "badge bg-danger",
        "Likely pathogenic": "badge bg-warning",
        "Uncertain": "badge bg-primary"
    }
    for idx, sig in df["Significance"].items():
        for significance, class_name in significance_badges.items():
            if sig.startswith(significance):
                df_for_table.loc[idx, "Significance"] = str(span(sig, cls=class_name))
                break
    return (df_for_table)
