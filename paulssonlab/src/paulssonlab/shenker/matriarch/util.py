import pandas as pd
from collections import defaultdict
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


def recursive_getattr(obj, keys):
    for k in keys:
        obj = obj[k]
    return obj


def diagnostics_to_dataframe(diagnostics):
    d = {
        k: flatten(v, predicate=lambda _, x: not isinstance(x, hv.ViewableElement))
        for k, v in diagnostics.items()
    }
    df = pd.DataFrame.from_dict(d, orient="index")
    return df


# FROM: https://stackoverflow.com/questions/6027558/flatten-nested-python-dictionaries-compressing-keys
def flatten(d, parent_key="", sep=".", predicate=None):
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + str(k) if parent_key else str(k)
        if isinstance(v, collections.MutableMapping):
            items.extend(flatten(v, new_key, sep=sep, predicate=predicate).items())
        else:
            if predicate is None or predicate(k, v):
                items.append((new_key, v))
    return dict(items)


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


def RevImage(img, **kwargs):
    return hv.Image(img[::-1], bounds=(0, 0, img.shape[1], img.shape[0])).opts(
        plot={"invert_yaxis": True}
    )


def RevRGB(img, **kwargs):
    return hv.RGB(img[::-1], bounds=(0, 0, img.shape[1], img.shape[0])).opts(
        plot={"invert_yaxis": True}
    )
