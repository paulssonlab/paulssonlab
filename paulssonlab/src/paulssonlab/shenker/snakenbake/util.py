import numpy as np

# from functools import lru_cache, wraps
import toolz
import cytoolz
from cytoolz import partial, compose


def make_odd(number):
    if number % 2 == 0:
        return number - 1
    else:
        return number


def make_hashable(obj):
    if isinstance(obj, np.ndarray):
        return hash(obj.tobytes())
    else:
        return toolz.sandbox.EqualityHashKey(None, obj)


def hash_key(args, kwargs):
    # return (args, hash(frozenset(kwargs.items())))
    # return (map(make_hashable, args), frozenset(kwargs.items()))
    args = tuple(map(make_hashable, args))
    kwargs = frozenset(map(compose(tuple, partial(map, make_hashable)), kwargs.items()))
    # print('args', args)
    # print('kwargs', kwargs)
    return (args, kwargs)
    # return (tuple(map(make_hashable, args)), frozenset(map(map(make_hashable), kwargs.items())))
    # return (map(make_hashable, args), frozenset(itemmap(map(make_hashable), kwargs)))


# memoize = lambda func: cachetools.cached({}, key=hashkey)(func)
memoize = lambda func: cytoolz.memoize(key=hash_key)(func)
