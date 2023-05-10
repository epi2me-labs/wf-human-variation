"""Common variables."""
from pandas.api import types as pd_types

CATEGORICAL = pd_types.CategoricalDtype(ordered=True)
CHROMOSOMES = {str(i): i for i in range(1, 23)}
CHROMOSOMES.update({f"chr{i}": int(i) for i in CHROMOSOMES})
CHROMOSOMES.update({'X': 23, 'Y': 24, 'chrX': 23, 'chrY': 24})
THEME = "epi2melabs"
