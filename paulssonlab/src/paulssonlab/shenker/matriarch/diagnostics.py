import pandas as pd
import holoviews as hv
from collections import defaultdict
from functools import wraps
from util import (
    wrap_dict_values,
    drop_dict_nones,
    flatten_dict,
    map_collections,
    tree,
    get_one,
)


def wrap_diagnostics(func, ignore_exceptions=False, pandas=False):
    @wraps(func)
    def wrapper(*args, **kwargs):
        diag = tree()
        # TODO: replace with util.fail_silently or toolz.excepts??
        err = None
        if ignore_exceptions:
            try:
                result = func(*args, **{"diagnostics": diag, **kwargs})
            except Exception as e:
                result = None
                err = e
        else:
            result = func(*args, **{"diagnostics": diag, **kwargs})
        if pandas:
            diag = diagnostics_to_series(diag)
        return (result, diag, err)

    return wrapper


def wrap_diagnostics_stack(func, **kwargs):
    func_diag = wrap_diagnostics(func)

    @wraps(func)
    def wrapper(img_stack, *args, **kwargs):
        return zip(*[func_diag(img, *args, **kwargs) for img in img_stack])

    return wrapper


def diagnostics_to_series(diagnostics):
    d = flatten_dict(
        diagnostics, predicate=lambda _, x: not isinstance(x, hv.ViewableElement)
    )
    df = pd.Series(d)
    return df


# TODO: replace with pd.melt??? probably not possible.
def expand_diagnostics_by_label(df, label="label_", keep_all=True):
    data = {}
    column_to_label_num = {}
    column_to_flattened_name = {}
    label_nums = set()
    for column in df.columns:
        column_parts = column.split(".")
        if column_parts[0].startswith(label):
            label_num = int(column_parts[0][len(label) :])
            column_to_label_num[column] = label_num
            column_to_flattened_name[column] = ".".join(column_parts[1:])
            label_nums.add(label_num)
    for idx, row in df.iterrows():
        d = defaultdict(dict)
        for column in df.columns:
            if column in column_to_label_num:
                label_num = column_to_label_num[column]
                d[label_num][column_to_flattened_name[column]] = row[column]
            else:
                if keep_all:
                    for label_num in label_nums:
                        d[label_num][column] = row[column]
        for label_num, dd in d.items():
            if isinstance(df.index, pd.MultiIndex):
                new_idx = idx + (label_num,)
            else:
                new_idx = (idx, label_num)
            data[new_idx] = dd
    df2 = pd.DataFrame.from_dict(data, orient="index")
    df2.index.rename(df.index.names + [label[:-1]], inplace=True)
    return df2


# TODO: unused??
def extract_diagnostics_singleton(diagnostics, keys):
    return wrap_dict_values(
        drop_dict_nones({k: get_in(keys, v) for k, v in diagnostics.items()}),
        ".".join(keys),
    )
