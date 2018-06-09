import pandas as pd
from collections import defaultdict, Sequence, Mapping
from tqdm import tqdm, tqdm_notebook
import zarr
from datetime import datetime, timezone
import holoviews as hv
from functools import reduce, partial, wraps
import wrapt
from cytoolz import compose
import collections
from dask.distributed import Future
import operator

# TODO: replace with toolz.excepts??
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


# TODO: replace with dicttoolz.get_in
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


# FROM: https://stackoverflow.com/questions/16458340/python-equivalent-of-zip-for-dictionaries
def zip_dicts(*dicts):
    for k in set(dicts[0]).intersection(*dicts[1:]):
        yield (k,) + tuple(d[k] for d in dicts)


# TODO: used?? made redundant by recursive_map??
def map_collections(func, data, max_level, contents=False):
    # TODO: allow specifying a range of levels
    if max_level == 0:
        if contents:
            return func(data)
        else:
            return data
    if isinstance(data, dict):
        d = {
            k: map_collections(func, v, max_level=max_level - 1)
            for k, v in data.items()
        }
        if contents:
            return d
        else:
            return func(d)
    else:
        raise NotImplementedError  # TODO: flesh out this function


# FROM: https://stackoverflow.com/questions/42095393/python-map-a-function-over-recursive-iterables/42095505
# TODO: document!!!
def recursive_map(func, data, shortcircuit=(), ignore=(), keys=False, max_level=None):
    if max_level is not None and max_level is not False:
        if max_level == 0:
            return func(data)
        max_level -= 1
    apply = lambda x: recursive_map(
        func,
        x,
        shortcircuit=shortcircuit,
        ignore=ignore,
        keys=keys,
        max_level=max_level,
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


map_futures = partial(recursive_map, shortcircuit=Future)

getitem_r = lambda b, a: operator.getitem(a, b)
getattr_r = lambda b, a: operator.getattr(a, b)


class Pointer:
    def __init__(self, value):
        self.value = value


def gather_futures(client, data):
    pointers = []
    futures = []

    def add_to_gather_list(future):
        pointer = Pointer(len(futures))
        pointers.append(pointer)
        futures.append(future)
        return pointer

    data_with_pointers = recursive_map(add_to_gather_list, data, shortcircuit=Future)
    results = client.gather(futures)
    pointer_to_result = dict(zip(pointers, results))
    return recursive_map(
        lambda p: pointer_to_result[p], data_with_pointers, shortcircuit=Pointer
    )


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
