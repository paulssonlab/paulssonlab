import numpy as np
import zarr
import re
import os
from cytoolz import partial, compose
from filelock import SoftFileLock
from contextlib import contextmanager
from collections import MutableMapping
from numcodecs import Blosc
import wrapt
from workflow import get_nd2_frame
from util import get_one, tqdm_auto, flatten_dict

DEFAULT_COMPRESSOR = Blosc(cname="zstd", clevel=5, shuffle=Blosc.SHUFFLE, blocksize=0)
DEFAULT_ORDER = "C"


def iterate_over_groupby(columns, progress_bar=tqdm_auto):
    @wrapt.decorator
    def wrapper(wrapped, instance, args, kwargs):
        res = {}
        pbar = progress_bar(args[0].groupby(columns))
        for key, group in pbar:
            key_kwargs = dict(zip(columns, key))
            pbar.set_postfix(key_kwargs)
            res[key] = wrapped(group, *args[1:], **{**kwargs, **key_kwargs})
        return res

    return wrapper


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
    # selected_trenches = trenches.xs((filename, position), drop_level=False)
    # TODO: here is where we specify logic for finding trenches for a given timepoint
    selected_trenches = trenches.xs((filename, position, 0), drop_level=False)
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
        for trench, trench_channels in channel_crops.items():
            for channel, crop in trench_channels.items():
                trench_crops.setdefault(trench, {}).setdefault(channel, {})[t] = crop
    for trench in trench_crops:
        if transformation is not None:
            trench_crops[trench] = {
                channel: transformation(crops)
                for channel, crops in trench_crops[trench].items()
            }
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


@contextmanager
def locked_zarr(
    filename, lock_dir=None, lock_suffix=".lock", zarr_store=zarr.LMDBStore
):
    # filename = '/tmp/foo.zarr'
    # zarr_store = zarr.ZipStore
    zarr_store = partial(zarr.LMDBStore, map_async=True)
    # zarr_store = zarr.DirectoryStore
    if lock_dir:
        lock_filename = os.path.join(lock_dir, os.path.basename(filename) + lock_suffix)
    else:
        lock_filename = filename + lock_suffix
    lock_dir = os.path.dirname(lock_filename)
    if not os.path.isdir(lock_dir):
        os.makedirs(lock_dir)
    lock = SoftFileLock(lock_filename)  # TODO: replace with fasteners??
    # import contextlib
    # lock = contextlib.nullcontext()
    with lock:
        store = zarr_store(filename)
        try:
            root = zarr.group(store=store, overwrite=False)
            yield root
        except:
            if hasattr(store, "close"):
                store.close()
            raise


def _zarr_filename(x, prefix=""):
    filename_part = (
        re.sub(".nd$", "", x["filename"], flags=re.IGNORECASE) + "." + prefix + ".zarr"
    )
    pos_part = "pos" + str(x["position"])
    return os.path.join(filename_part, pos_part)


DEFAULT_OUTPUT_FUNC = compose(locked_zarr, _zarr_filename)

_values_are_not_dict = lambda x: isinstance(get_one(x), MutableMapping)


def write_images_to_zarr(
    image_data,
    root,
    merge=True,
    overwrite=False,
    compressor=DEFAULT_COMPRESSOR,
    order=DEFAULT_ORDER,
):
    for key, t_images in flatten_dict(
        image_data, lookahead=_values_are_not_dict
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


# overwrite=True: write to temp path, atomic move into position after write
# merge=True: read in, slot in
# overwrite=False: raise error if zarr exists
@iterate_over_groupby(["filename", "position"])
def compress_frames(
    frames,
    trenches,
    output_func=DEFAULT_OUTPUT_FUNC,
    merge=True,
    overwrite=False,
    compressor=DEFAULT_COMPRESSOR,
    order=DEFAULT_ORDER,
    filename=None,
    position=None,
    **kwargs,
):
    trench_crops = _get_trench_crops(
        frames, trenches, filename=filename, position=position, **kwargs
    )
    with output_func(dict(filename=filename, position=position)) as root:
        write_images_to_zarr(
            trench_crops,
            root,
            merge=merge,
            overwrite=overwrite,
            compressor=compressor,
            order=order,
        )
    return trench_crops
