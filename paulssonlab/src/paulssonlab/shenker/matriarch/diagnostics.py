import pandas as pd
import holoviews as hv
from collections import defaultdict
from functools import wraps
from util import wrap_dict_values, drop_dict_nones, flatten_dict, map_collections, tree


def wrap_diagnostics(func, ignore_exceptions=True):
    @wraps(func)
    def wrapper(*args, **kwargs):
        diag = tree()
        # TODO: replace with util.fail_silently or toolz.excepts??
        if ignore_exceptions:
            try:
                result = func(*args, **{"diagnostics": diag, **kwargs})
            except Exception as e:
                return (None, diag, e)
        else:
            result = func(*args, **{"diagnostics": diag, **kwargs})
        return (result, diag, None)

    return wrapper


# @wrapt.decorator
# def wrap_diagnostics(func, instance, args, kwargs):
#     diag = tree()
#     return (func(*args, **{'diagnostics': diag, **kwargs}), diag)


def diagnostics_to_dataframe(diagnostics):
    d = flatten_dict(
        diagnostics, predicate=lambda _, x: not isinstance(x, hv.ViewableElement)
    )
    df = pd.DataFrame.from_dict({0: d}, orient="index")
    return df


def wrapped_diagnostics_to_dataframe(x):
    return expand_diagnostics_by_label(diagnostics_to_dataframe(x[1]))


# TODO: replace with pd.melt???
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
        # data.extend(d.values())
        for label_num, dd in d.items():
            data[(idx, label_num)] = dd
    df2 = pd.DataFrame.from_dict(data, orient="index")
    # parent_index = pd.MultiIndex(df.index)
    # df2.index.rename(parent_index.names, level=parent_index.levels, inplace=True)
    df2.index.rename([df.index.name, label[:-1]], inplace=True)
    return df2


def extract_diagnostics_singleton(diagnostics, keys):
    return wrap_dict_values(
        drop_dict_nones({k: get_in(keys, v) for k, v in diagnostics.items()}),
        ".".join(keys),
    )
