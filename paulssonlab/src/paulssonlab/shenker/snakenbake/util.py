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
        # mutable ndarrays aren't hashable
        return hash(obj.tobytes())
    elif isinstance(obj, list):
        # if we don't handle list, we get a weird error (TODO)
        return hash(tuple(obj))
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


# memoize = lambda func: cachetools.cached({}, key=hash_key)(func)
memoize = lambda func: cytoolz.memoize(key=hash_key)(func)
# memoize = lambda x: x
