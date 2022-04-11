import numpy as np
import pandas as pd


def valid_range(ary):
    is_valid = ~np.isnan(ary)
    idx1 = is_valid.argmax()
    idx2 = len(ary) - is_valid[::-1].argmax()
    return idx1, idx2


def bin_size(x):
    return (x[-1] - x[0]) / len(x)


def interpolate_dataframe(df, new_index, union=False):
    old_index = df.index.astype(np.float)
    if not isinstance(new_index, pd.core.indexes.base.Index):
        new_index = pd.Index(new_index, name=df.index.name)
    union_index = new_index.union(old_index)
    new_df = df.reindex(index=union_index)
    new_df.interpolate(method="linear", limit_area="inside", inplace=True)
    if not union:
        new_df = new_df.loc[new_index]
    return new_df


def bin_dataframe(df, bins):
    bin_assignment = pd.cut(df.index, bins).rename_categories(
        (bins.right + bins.left) / 2
    )
    return df.groupby(bin_assignment).mean()
