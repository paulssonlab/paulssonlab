# import wrapt # TODO
# from decorator import decorator
import functools
import gc
import os
import re
from contextlib import contextmanager

import numpy as np
import zarr
from cytoolz import compose, partial
from filelock import SoftFileLock
from numcodecs import Blosc
from tqdm.auto import tqdm

from paulssonlab.image_analysis.data_io import (
    write_dataframe_to_arrow,
    write_dataframe_to_parquet,
)
from paulssonlab.image_analysis.util import (
    flatten_dict,
    get_one,
    mapping_values_are_dict,
    unflatten_dict,
)
from paulssonlab.image_analysis.workflow import get_nd2_frame

DEFAULT_COMPRESSOR = Blosc(cname="zstd", clevel=5, shuffle=Blosc.SHUFFLE, blocksize=0)
DEFAULT_ORDER = "C"

DEFAULT_LOCK_TIMEOUT = 120  # seconds

zarrify = partial(
    zarr.array, compressor=DEFAULT_COMPRESSOR, order=DEFAULT_ORDER, chunks=False
)

# @functools.decorator
# def spawn(func, *args, **kwargs):
#     # TODO: use queues to return results
#     t = multiprocessing.Process(target=func, args=args, kwargs=kwargs)
#     t.start()
#     t.join()

# def iterate_over_groupby(columns, progress_bar=tqdm):
#     # TODO
#     @wrapt.decorator
#     def wrapper(wrapped, instance, args, kwargs):
#     #@decorator
#     #@functools.wraps
#     #def wrapper(wrapped, *args, **kwargs):
#         res = {}
#         pbar = progress_bar(args[0].groupby(columns))
#         for key, group in pbar:
#             key_kwargs = dict(zip(columns, key))
#             pbar.set_postfix(key_kwargs)
#             res[key] = wrapped(group, *args[1:], **{**kwargs, **key_kwargs})
#         return res
#     return wrapper


def iterate_over_groupby(columns, progress_bar=None):  # tqdm):
    def get_wrapper(wrapped):
        @functools.wraps(wrapped)
        def wrapper(*args, **kwargs):
            res = {}
            groups = args[0].groupby(columns)
            if progress_bar is not None:
                groups = progress_bar(groups)
            for key, group in groups:
                key_kwargs = dict(zip(columns, key))
                if progress_bar is not None:
                    groups.set_postfix(key_kwargs)
                res[key] = wrapped(group, *args[1:], **{**kwargs, **key_kwargs})
            return res

        # wrapper.__name__ == wrapped.__name__
        # wrapper.__doc__ == wrapped.__doc__
        # wrapper.__module__ == wrapped.__module__
        return wrapper

    return get_wrapper


def index_dict_to_array(d, array_func=np.zeros):
    idx_max = max(d.keys())
    shape0 = get_one(d).shape
    shape = shape0 + (idx_max + 1,)
    ary = array_func(shape)
    for k, v in d.items():
        ary[..., k] = v
    return ary


def _get_trench_crops(
    trenches,
    frames,
    get_frame_func=get_nd2_frame,
    transformation=None,
    include_frame=False,
    frame_transformation=None,
    filename=None,
    position=None,
):
    # selected_trenches = trenches.groupby(['filename', 'position']).get_group((filename, position)
    selected_trenches = trenches.xs((filename, position), drop_level=False)
    # TODO: here is where we specify logic for finding trenches for a given timepoint
    # selected_trenches = trenches.xs((filename, position, 0), drop_level=False)
    trench_crops = {}
    if include_frame:
        whole_frames = {}
    for t, channel_frames in frames.groupby("t"):
        channels = set(channel_frames.index.get_level_values("channel"))
        channel_images = {
            channel: get_frame_func(
                filename=filename, position=position, channel=channel, t=t
            )
            for channel in channels
        }
        if include_frame:
            for channel, image in channel_images.items():
                if frame_transformation is not None:
                    image = frame_transformation(image)
                whole_frames.setdefault(channel, {})[t] = image
        channel_crops = _crop_trenches(selected_trenches, channel_images)
        for trench_key, trench_channels in channel_crops.items():
            for channel, crop in trench_channels.items():
                trench_crops.setdefault(trench_key, {}).setdefault(channel, {})[
                    t
                ] = crop
        del channel_images
        gc.collect()
    for trench_key in trench_crops:
        if transformation is not None:
            trench_crops[trench_key] = {
                channel: transformation(crops)
                for channel, crops in trench_crops[trench_key].items()
            }
    # break out
    trench_crops = unflatten_dict(trench_crops)
    if include_frame:
        for channel, frames in whole_frames.items():
            trench_crops.setdefault("_frame", {})[channel] = frames
    return trench_crops


get_trench_crops = iterate_over_groupby(["filename", "position"])(_get_trench_crops)
_get_trench_stacks = partial(_get_trench_crops, transformation=index_dict_to_array)
get_trench_stacks = iterate_over_groupby(["filename", "position"])(_get_trench_stacks)

# def crop_trenches(trenches, image, filename=None,
#                   position=None, t=None):
#     #selected_trenches = trenches.groupby(['filename', 'position', 't']).get_group((filename, position, t)
#     # TODO: for now, only select trench definitions for first timepoint
#     selected_trenches = get_one(trenches.groupby(['filename', 'position']).get_group((filename, position)).groupby('t'))[1]
#     return _crop_trenches(selected_trenches, image)


def _crop_trenches(trenches, images):
    trench_sets = trenches.index.get_level_values("trench_set")
    trench_idxs = trenches.index.get_level_values("trench")
    # this requires re-allocating arrays unless columns were stored next to each other by pandas blockmanager
    # uls = trenches['upper_left'].values
    # lrs = trenches['lower_right'].values
    uls_x = trenches[("upper_left", "x")].values
    uls_y = trenches[("upper_left", "y")].values
    lrs_x = trenches[("lower_right", "x")].values
    lrs_y = trenches[("lower_right", "y")].values
    trench_crops = {}
    for i in range(len(uls_x)):
        ul_x = uls_x[i]
        ul_y = uls_y[i]
        lr_x = lrs_x[i]
        lr_y = lrs_y[i]
        channel_crops = {
            channel: channel_image[ul_y : lr_y + 1, ul_x : lr_x + 1]
            for channel, channel_image in images.items()
        }
        trench_crops[(trench_sets[i], trench_idxs[i])] = channel_crops
    return trench_crops


def file_lock(filename, lock_dir=None, lock_suffix=".lock", makedirs=True):
    if lock_dir:
        lock_filename = os.path.join(lock_dir, os.path.basename(filename) + lock_suffix)
    else:
        lock_filename = filename + lock_suffix
    if makedirs:
        lock_dir = os.path.dirname(lock_filename)
        os.makedirs(lock_dir, exist_ok=True)
        os.makedirs(os.path.dirname(filename), exist_ok=True)
    lock = SoftFileLock(lock_filename)  # TODO: replace with fasteners??
    return lock


@contextmanager
def open_zarr(filename, store=partial(zarr.LMDBStore, lock=False)):
    store_ = store(filename)
    try:
        root = zarr.group(store=store_, overwrite=False)
        yield root
    except:
        if hasattr(store_, "close"):
            store_.close()
        raise


def write_images_to_zarr(
    image_data,
    root,
    merge=True,
    overwrite=False,
    compressor=DEFAULT_COMPRESSOR,
    order=DEFAULT_ORDER,
):
    for key, t_images in flatten_dict(
        image_data, lookahead=mapping_values_are_dict
    ).items():
        arr_path = "/".join(map(str, key))
        dtype = get_one(t_images).dtype
        shape = get_one(t_images).shape
        if merge:
            old_arr = root.get(arr_path, None)
        else:
            old_arr = None
        t_max = max(t_images.keys())
        if old_arr is not None:
            t_max = max(old_arr.shape[-1] - 1, t_max)
        shape = shape + (t_max + 1,)
        chunks = False  # (1,1,1) # TODO: different for _frame
        # TODO: copy here is unnecessary unless merging
        # TODO: handle dtype
        new_data = np.zeros(shape)
        if old_arr is not None:
            new_data[:, :, : old_arr.shape[-1]] = old_arr
        for t, image in t_images.items():
            new_data[:, :, t] = image
        # print(f'writing array {filename} pos{position} ch:{channel} {trench_set}.{trench}.t{t}')
        if merge or overwrite:
            try:
                del root[arr_path]
            except KeyError:
                pass
        new_arr = root.create_dataset(
            arr_path,
            overwrite=False,
            data=new_data,
            shape=shape,
            dtype=dtype,
            compressor=compressor,
            chunks=chunks,
            order=order,
        )
    return image_data


def write_images_and_measurements(
    d,
    filename_func,
    merge=True,
    overwrite=True,
    write_images=True,
    write_measurements=True,
    dataframe_format="arrow",
    compressor=DEFAULT_COMPRESSOR,
    order=DEFAULT_ORDER,
    timeout=DEFAULT_LOCK_TIMEOUT,
):
    for key, results in d.items():
        key_dict = dict(zip(("filename", "position"), key))
        # write images
        if "measurements" in results and write_measurements:
            for name, df in results["measurements"].items():
                if dataframe_format == "arrow":
                    measurements_filename = filename_func(
                        kind="measurements", extension="arrow", name=name, **key_dict
                    )
                    write_func = write_dataframe_to_arrow
                elif dataframe_format == "parquet":
                    measurements_filename = filename_func(
                        kind="measurements", extension="parquet", name=name, **key_dict
                    )
                    write_func = write_dataframe_to_parquet
                else:
                    raise ValueError
                with file_lock(measurements_filename).acquire(timeout=timeout) as lock:
                    write_func(
                        measurements_filename, df, merge=merge, overwrite=overwrite
                    )
        # write images
        if "images" in results and write_images:
            images_filename = filename_func(kind="images", extension="zarr", **key_dict)
            with file_lock(images_filename).acquire(timeout=timeout) as lock, open_zarr(
                images_filename
            ) as root:
                write_images_to_zarr(
                    results["images"],
                    root,
                    merge=merge,
                    overwrite=overwrite,
                    compressor=compressor,
                    order=order,
                )
    return d


# overwrite=True: write to temp path, atomic move into position after write
# merge=True: read in, slot in
# overwrite=False: raise error if zarr exists
# @iterate_over_groupby(['filename', 'position'])
# def compress_frames(frames, trenches, output_func=DEFAULT_OUTPUT_FUNC,
#                     merge=True,
#                     overwrite=False,
#                     compressor=DEFAULT_COMPRESSOR,
#                     order=DEFAULT_ORDER,
#                     filename=None,
#                     position=None,
#                     **kwargs):
#     trench_crops = _get_trench_crops(frames, trenches,
#                                      filename=filename,
#                                      position=position,
#                                      **kwargs)
#     with output_func(dict(filename=filename, position=position)) as root:
#         write_images_to_zarr(trench_crops,
#                              root,
#                              merge=merge,
#                              overwrite=overwrite,
#                              compressor=compressor,
#                              order=order)
#     return trench_crops
