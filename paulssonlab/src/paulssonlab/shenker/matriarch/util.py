from collections import defaultdict
from tqdm import tqdm, tqdm_notebook
import zarr
from datetime import datetime, timezone
import holoviews as hv


def fail_silently(func):
    try:
        return func()
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
    return tqdm(*args, **kwargs)
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
