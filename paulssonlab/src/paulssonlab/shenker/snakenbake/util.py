import numpy as np

# from functools import lru_cache, wraps
import toolz
import cytoolz
from cytoolz import partial, compose
import gdspy
import matplotlib.pyplot as plt


def plot_cell(cell, exclude=(2,)):
    # FROM: https://github.com/heitzmann/gdspy/issues/42
    if hasattr(cell, "get_polygons"):
        poly_dict = cell.get_polygons(by_spec=True)
    else:
        poly_dict = {
            l_d: [p] for l_d, p in zip(zip(cell.layers, cell.datatypes), cell.polygons)
        }
    plt.figure(figsize=(40, 20))
    for layer_datatype, polys in poly_dict.items():
        if layer_datatype[0] in exclude:
            continue
        for poly in polys:
            plt.fill(*poly.T, lw=0.5, ec="k", fc=(1, 0, 0, 0.5))
    plt.axes().set_aspect("equal", "datalim")


def write_gds(main_cell, filename, unit=1.0e-6, precision=1.0e-9):
    cells = [main_cell] + list(main_cell.get_dependencies(True))
    gdspy.write_gds(filename, cells=cells, unit=unit, precision=precision)


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
