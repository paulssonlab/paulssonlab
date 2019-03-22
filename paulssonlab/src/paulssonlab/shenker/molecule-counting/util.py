import numpy as np


def short_circuit_none(func, *args, **kwargs):
    if args[-1] is None:
        return None
    else:
        return func(*args, **kwargs)


# TODO: unused
def trim_nones(ary):
    length = len(ary)
    for idx, elem in enumerate(reversed(ary)):
        if elem is not None:
            length = idx
            break
    return ary[: len(ary) - length]


def none_to_nans(iterable):
    iterable = list(iterable)
    non_none = None
    for elem in iterable:
        if elem is not None:
            non_none = elem
            break
    if non_none is None:
        return None
    replacement = np.full(np.shape(non_none), np.nan)
    return [elem if elem is not None else replacement for elem in iterable]
