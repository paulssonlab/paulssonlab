from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass
from functools import partial
from itertools import zip_longest
from typing import Callable

import dask
from cytoolz import compose, reduce, take
from distributed import Client


def tree():
    return defaultdict(tree)


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


def getitem_if_not_none(obj, key):
    if obj is not None:
        return obj[key]
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


@dataclass
class DaskDistributedDelayed:
    func: Callable
    client: Client

    def __call__(self, *args, **kwargs):
        return self.client.submit(self.func, *args, **kwargs)


def get_delayed(delayed):
    if delayed is True:
        return dask.delayed(pure=True)
    elif delayed is False:

        def wrapper(func, **kwargs):
            return func

        return wrapper
    elif isinstance(delayed, Client):

        def outer_wrapper(func, **kwargs):
            return DaskDistributedDelayed(func, delayed)

        return outer_wrapper
    else:
        return delayed


def repeat_apply(func, n):
    if n <= 0:
        return lambda x: x
    return reduce(lambda f1, f2: compose(f1, f2), [func] * n)
