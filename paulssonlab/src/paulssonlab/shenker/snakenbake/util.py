import numpy as np

# from functools import lru_cache, wraps
import toolz
import cytoolz
from cytoolz import partial, compose
import gdstk
import matplotlib.pyplot as plt
import shortuuid

get_uuid = partial(shortuuid.random, length=4)


def plot_cell(polys, exclude=(2,)):
    plt.figure(figsize=(10, 5), dpi=300)
    plt.axes().set_aspect("equal", "datalim")
    if hasattr(polys, "get_polygons"):
        polys = polys.get_polygons()
    elif hasattr(polys, "points"):
        polys = [polys]
    else:
        raise NotImplementedError
    for poly in polys:
        plt.fill(*poly.points.T, lw=0.5, ec="k", fc=(1, 0, 0, 0.5))


def write_gds(main_cell, filename, unit=1.0e-6, precision=1.0e-9, max_points=199):
    cells = [main_cell] + list(main_cell.dependencies(True))
    lib = gdstk.Library(unit=unit, precision=precision)
    lib.add(*cells)
    lib.write_gds(filename, max_points)


def strip_units(x):
    if hasattr(x, "magnitude"):
        return x.magnitude
    else:
        return x


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
