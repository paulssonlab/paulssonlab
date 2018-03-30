import pandas as pd
from collections import defaultdict, Sequence, Mapping
from tqdm import tqdm, tqdm_notebook
import zarr
from datetime import datetime, timezone
import holoviews as hv
from functools import reduce
from cytoolz import compose
import collections


def fail_silently(func):
    try:
        return func()
    except KeyboardInterrupt:
        raise
    except:
        return None


def getattr_if_not_none(obj, key):
    if obj is not None:
        return obj[key]
    else:
        return None


def iterate_getattr(obj, keys):
    for k in keys:
        if k in obj:
            obj = obj[k]
        else:
            return None
    return obj


def iterate_get_collection_value(obj, level):
    for i in range(level):
        if isinstance(obj, Mapping):
            obj = iter(obj.values())
        obj = next(obj)
    return obj


def extract_diagnostics_singleton(diagnostics, keys):
    return wrap_dict_values(
        drop_dict_nones({k: iterate_getattr(v, keys) for k, v in diagnostics.items()}),
        ".".join(keys),
    )


def wrap_dict_values(obj, label):
    return {k: {label: v} for k, v in obj.items()}


def drop_dict_nones(obj):
    return {k: v for k, v in obj.items() if v is not None}


# FROM: https://stackoverflow.com/questions/4664850/find-all-occurrences-of-a-substring-in-python
def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1:
            return
        yield start
        start += 1  # find overlapping matches (otherwise += len(sub))


def diagnostics_to_dataframe(diagnostics):
    d = {
        k: flatten_dict(v, predicate=lambda _, x: not isinstance(x, hv.ViewableElement))
        for k, v in diagnostics.items()
    }
    df = pd.DataFrame.from_dict(d, orient="index")
    df.index.rename("pos", inplace=True)
    return df


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


def drop_constant_columns(df):
    return df.loc[:, df.nunique(dropna=False) != 1]


# FROM: https://stackoverflow.com/questions/6027558/flatten-nested-python-dictionaries-compressing-keys
def flatten_dict(d, parent_key="", sep=".", predicate=None):
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + str(k) if parent_key else str(k)
        if isinstance(v, collections.MutableMapping):
            items.extend(flatten_dict(v, new_key, sep=sep, predicate=predicate).items())
        else:
            if predicate is None or predicate(k, v):
                items.append((new_key, v))
    return dict(items)


def map_collections(func, data, max_level):
    # TODO: allow specifying a range of levels
    if max_level == 0:
        return data
    if isinstance(data, dict):
        return func(
            {
                k: map_collections(func, v, max_level=max_level - 1)
                for k, v in data.items()
            }
        )
    else:
        raise NotImplementedError  # TODO: flesh out this function


# FROM: https://stackoverflow.com/questions/42095393/python-map-a-function-over-recursive-iterables/42095505
# TODO: document!!!
def recursive_map(func, data, shortcircuit=(), ignore=(), keys=False):
    apply = lambda x: recursive_map(
        func, x, shortcircuit=shortcircuit, ignore=ignore, keys=keys
    )
    if isinstance(data, shortcircuit):
        return func(data)
    # to avoid an infinite recursion error
    # we need to hardcode that strs are ignored, because str[0] is a str, and hence a Sequence
    elif isinstance(data, str) or ignore is not True and isinstance(data, ignore):
        return data
    elif isinstance(data, Mapping):
        if keys:
            return type(data)(dict({apply(k): apply(v) for k, v in data.items()}))
        else:
            return type(data)(dict({k: apply(v) for k, v in data.items()}))
    elif isinstance(data, Sequence):
        return type(data)(list(apply(v) for v in data))
    elif ignore is True or True in ignore:
        return data
    else:
        return func(data)


def repeat_apply(func, n):
    if n <= 0:
        return lambda x: x
    return reduce(lambda f1, f2: compose(f1, f2), [func] * n)


# FROM: https://stackoverflow.com/questions/2150739/iso-time-iso-8601-in-python
def isonow():
    return datetime.datetime.utcnow().replace(tzinfo=datetime.timezone.utc).isoformat()


# FROM: https://stackoverflow.com/questions/2150739/iso-time-iso-8601-in-python
def timestamp_to_isoformat(ts):
    dt = datetime.fromtimestamp(ts, timezone.utc)
    return dt.astimezone().isoformat()


def tree():
    return defaultdict(tree)


def tqdm_auto(*args, **kwargs):
    try:
        # FROM: https://github.com/tqdm/tqdm/issues/372
        from IPython import get_ipython

        ipy = get_ipython()
        if ipy and ipy.__class__.__name__ == "ZMQInteractiveShell":
            return tqdm_notebook(*args, **kwargs)
    except:
        pass
    return tqdm(*args, **kwargs)


# return tqdm(*args, **kwargs)
# try:
#    return tqdm_notebook(*args, **kwargs)
# except:
#    return tqdm(*args, **kwargs)


def open_zarr_group(dir_path):
    store = zarr.DirectoryStore(dir_path)
    return zarr.open_group(store=store)
