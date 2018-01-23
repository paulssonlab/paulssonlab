from collections import defaultdict
from tqdm import tqdm, tqdm_notebook
import zarr


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
