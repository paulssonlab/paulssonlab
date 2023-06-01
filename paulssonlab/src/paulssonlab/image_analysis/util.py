import collections
import functools
import numbers
import operator
import os
import random
from collections import defaultdict, namedtuple
from collections.abc import Mapping, MutableMapping, Sequence
from datetime import datetime, timezone
from functools import wraps
from itertools import count, zip_longest

import holoviews as hv
import numpy as np
import pandas as pd
import pandas.core.groupby.groupby
import wrapt
import zarr
from cytoolz import compose, keymap, partial, reduce, take
from dask.distributed import Future
from tqdm import tqdm, tqdm_notebook


def return_none(*args, **kwargs):
    return None


# FROM: http://louistiao.me/posts/adding-__name__-and-__doc__-attributes-to-functoolspartial-objects/
def wrapped_partial(func, *args, **kwargs):
    partial_func = functools.partial(func, *args, **kwargs)
    functools.update_wrapper(partial_func, func)
    return partial_func


def format_number(x, digits):
    if isinstance(x, numbers.Integral):
        return str(x)
    else:
        return "{0:.{1}f}".format(value, digits)


def kwcompose(func1, func2):
    @wraps(func2)
    def func(*args, **kwargs):
        return func1(**func2(*args, **kwargs))

    return func


# FROM: https://stackoverflow.com/questions/23937433/efficiently-joining-two-dataframes-based-on-multiple-levels-of-a-multiindex?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
def multi_join(left, right, **kwargs):
    if isinstance(left, pd.Index):
        left_mergable = left.to_frame(index=False)
        left_index = left
    elif isinstance(left, pd.DataFrame):
        left_mergable = left.reset_index()
        left_index = left.index
    else:
        raise NotImplementedError
    res = pd.merge(
        left_mergable,
        right.reset_index(),
        **{"on": right.index.names, "how": "inner", **kwargs},
    ).set_index(left_index.names)
    if isinstance(right, pd.Series):
        res = res.squeeze()
    return res


# def iter_index(df):
#     if isinstance(df, pandas.core.groupby.groupby.GroupBy):
#         Index = namedtuple('Index', df.keys, rename=True) # TODO: have not tested rename=True
#         for idx, group in df:
#             yield Index(*idx), group
#     else:
#         if isinstance(df, pd.Index):
#             index = df
#             rows = df.to_frame(index=False)
#         else:
#             index = df.index
#             rows = df.reset_index()
#         # the following line will be 4-6x faster in Python 3.7
#         # CF: https://bugs.python.org/issue28638
#         Index = namedtuple('Index', index.names, rename=True) # TODO: have not tested rename=True
#         for row in rows.itertuples():
#             idx = Index(*index[row.Index])
#             yield idx, row


def split_into(ary, max_length):
    if max_length is None:
        return ary
    return np.array_split(ary, max(len(ary) // max_length, 1))


UNEVEN_GROUPS = object()


# FROM: https://docs.python.org/3/library/itertools.html#recipes
def grouper(iterable, n, fillvalue=UNEVEN_GROUPS):
    """Collect data into fixed-length chunks or blocks."""
    args = [iter(iterable)] * n
    filter_filler = False
    if fillvalue == UNEVEN_GROUPS:
        fillvalue = object()
        filter_filler = True
    groups = zip_longest(*args, fillvalue=fillvalue)
    if filter_filler:
        groups = map(compose(tuple, partial(filter, lambda x: x != fillvalue)), groups)
    return groups


# TODO: replace with toolz.excepts??
def fail_silently(func):
    try:
        return func()
    except KeyboardInterrupt:
        raise
    except:
        return None


def getitem_if_not_none(obj, key):
    if obj is not None:
        return obj[key]
    else:
        return None


def getattr_if_not_none(obj, key):
    if obj is not None:
        return getattr(obj, key)
    else:
        return None


def _get_one(obj, count=1):
    if isinstance(obj, Mapping):
        obj = obj.values()
    res = list(take(count, iter(obj)))
    if count == 1:
        return res[0]
    else:
        return res


def get_one(obj, count=1, level=1):
    return _get_one(repeat_apply(_get_one, level - 1)(obj), count=count)


get_some = partial(get_one, count=3)


def get_random(obj, count=1, level=1):
    if isinstance(obj, Mapping):
        obj = obj.values()
    obj = list(obj)
    obj = repeat_apply(get_random, level - 1)(obj)
    res = random.sample(obj, count)
    if count == 1:
        return res[0]
    else:
        return res


def get_keys(d, *__keys, keys=None):
    if keys:
        keys = set(keys)
        if __keys:
            raise ValueError("only one of __keys or keys can be specified")
    else:
        keys = set(__keys)
    return {k: v for k, v in d.items() if k in keys}


def get_kwargs(*__keys, keys=None, **kwargs):
    return get_keys(kwargs, *__keys, keys=keys)


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


# predicate that's useful for flatten_dict
mapping_values_are_dict = lambda x: isinstance(get_one(x), MutableMapping)


# FROM: https://stackoverflow.com/questions/6027558/flatten-nested-python-dictionaries-compressing-keys
def flatten_dict(d, parent_key=None, sep=None, predicate=None, lookahead=None):
    if parent_key is None:
        if sep is not None:
            parent_key = ""
        else:
            parent_key = ()
    items = []
    for k, v in d.items():
        if sep is not None:
            new_key = parent_key + sep + str(k) if parent_key else str(k)
        else:
            new_key = parent_key + (k,)
        if isinstance(v, MutableMapping) and (lookahead is None or lookahead(v)):
            items.extend(
                flatten_dict(
                    v,
                    parent_key=new_key,
                    sep=sep,
                    predicate=predicate,
                    lookahead=lookahead,
                ).items()
            )
        else:
            if predicate is None or predicate(k, v):
                items.append((new_key, v))
    return dict(items)


# FROM: https://stackoverflow.com/questions/50607128/creating-a-nested-dictionary-from-a-flattened-dictionary
def unflatten_dict(d, sep=None):
    output = {}
    for k, v in d.items():
        current = output
        if sep is not None:
            pieces = k.split(sep)
        else:
            pieces = k
        for piece in pieces[:-1]:
            if not piece in current:
                current[piece] = {}
            current = current[piece]
        current[pieces[-1]] = v
    return output


def unzip_items(d):
    return list(zip(*list(d)))


# FROM: https://stackoverflow.com/questions/16458340/python-equivalent-of-zip-for-dictionaries
def zip_dicts(*dicts):
    for k in set(dicts[0]).intersection(*dicts[1:]):
        yield (k,) + tuple(d[k] for d in dicts)


# simple one-level version of unzip_collection
def unzip_dicts(d):
    val0 = get_one(d.values())
    length = len(val0)
    return [{k: v[idx] for k, v in d.items()} for idx in range(length)]


def map_dict_levels(func, d):
    return unflatten_dict(keymap(compose(tuple, func), flatten_dict(d)))


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
def recursive_map(
    func,
    data,
    shortcircuit=(),
    ignore=(),
    keys=False,
    max_level=None,
    predicate=None,
    key_predicate=None,
):
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
        predicate=predicate,
        key_predicate=key_predicate,
    )
    if isinstance(data, shortcircuit):
        return func(data)
    # to avoid an infinite recursion error
    # we need to hardcode that strs are ignored, because str[0] is a str, and hence a Sequence
    elif isinstance(data, str) or ignore is not True and isinstance(data, ignore):
        return data
    elif isinstance(data, Mapping):
        if keys:
            values = {apply(k): apply(v) for k, v in data.items()}
        else:
            values = {k: apply(v) for k, v in data.items()}
        if key_predicate is not None:
            values = {k: v for k, v in values.items() if key_predicate(k)}
        if predicate is not None:
            values = {k: v for k, v in values.items() if predicate(v)}
        return type(data)(values)
    elif isinstance(data, Sequence):
        values = [apply(v) for v in data]
        if predicate is not None:
            values = [v for v in values if predicate(v)]
        return type(data)(values)
    elif ignore is True or True in ignore:
        return data
    else:
        return func(data)


map_futures = partial(recursive_map, shortcircuit=Future)

getitem_r = lambda b, a: operator.getitem(a, b)
getattr_r = lambda b, a: operator.getattr(a, b)


# TODO: unused. remove? see unzip_dicts, above
def unzip_collection(data, **kwargs):
    for idx in count():
        yield recursive_map(
            partial(getitem_r, idx), data, **{"shortcircuit": tuple, **kwargs}
        )


class Pointer:
    def __init__(self, value):
        self.value = value


def apply_map_futures(func, data, predicate=None, skip_empty=True):
    pointers = []
    futures = []
    pointer_num = count()
    skip_value = object()

    def add_to_future_list(future):
        pointer = Pointer(next(pointer_num))
        if predicate is None or predicate(future):
            pointers.append(pointer)
            futures.append(future)
        return pointer

    data_with_pointers = map_futures(add_to_future_list, data)
    results = func(futures)
    pointer_to_result = dict(zip(pointers, results))
    if predicate:
        if skip_empty:
            skip_predicate = lambda x: x is not skip_value and (
                not hasattr(x, "__len__") or len(x) > 0
            )
        else:
            skip_predicate = lambda x: x is not skip_value
    else:
        skip_predicate = None
    return recursive_map(
        lambda p: pointer_to_result.get(p, skip_value),
        data_with_pointers,
        shortcircuit=Pointer,
        predicate=skip_predicate,
    )


# TODO: rewrite all _future util functions in terms of generic shortcutting recursive_map functions
def collect_futures(data, predicate=lambda f: True):
    futures = []
    map_futures(lambda f: futures.append(f) if predicate(f) else None, data)
    return futures


finished_futures = partial(collect_futures, predicate=lambda f: f.status == "finished")
pending_futures = partial(collect_futures, predicate=lambda f: f.status == "pending")
failed_futures = partial(collect_futures, predicate=lambda f: f.status == "error")


def array_to_tuples(ary):
    return [tuple(a) for a in ary]


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


def open_zarr_group(dir_path):
    store = zarr.DirectoryStore(dir_path)
    return zarr.open_group(store=store)


def summarize_filenames(filenames):
    common_prefix = os.path.commonpath(filenames)
    return ["/..." + filename[len(common_prefix) :] for filename in filenames]
