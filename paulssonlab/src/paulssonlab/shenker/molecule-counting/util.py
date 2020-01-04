import numpy as np
import dask.base


def conditional(flag, true_val, false_val):
    if flag:
        return true_val
    else:
        return false_val


# TODO: unused
def short_circuit_first(func, *args, **kwargs):
    if args[0] is None:
        return None
    else:
        return func(*args, **kwargs)


# TODO: unused
def short_circuit_last(func, *args, **kwargs):
    if args[-1] is None:
        return None
    else:
        return func(*args, **kwargs)
