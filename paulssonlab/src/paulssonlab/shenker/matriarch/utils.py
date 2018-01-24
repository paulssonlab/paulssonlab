from collections import defaultdict
from tqdm import tqdm, tqdm_notebook
import zarr
from datetime import datetime, timezone

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
