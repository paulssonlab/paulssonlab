from itertools import product
from numbers import Integral

import cachetools
import dask
import h5py
import nd2reader
import numpy as np
from cytoolz import excepts
from tqdm.auto import tqdm

from paulssonlab.io.metadata import parse_nd2_metadata

ND2READER_CACHE = cachetools.LFUCache(maxsize=48)


def _get_nd2_reader(filename, **kwargs):
    # TODO: removed memmap flag
    return nd2reader.ND2Reader(filename, **kwargs)


get_nd2_reader = cachetools.cached(cache=ND2READER_CACHE)(_get_nd2_reader)


def get_nd2_frame(filename, position, channel_num, t):
    reader = get_nd2_reader(filename)
    ary = reader.get_frame_2D(v=position, c=channel_num, t=t)
    return ary


def _select_indices(idxs, slice_):
    if slice_ is None:
        return idxs
    elif isinstance(slice_, slice):
        return idxs[slice_]
    elif isinstance(slice_, Integral):
        return idxs[:slice_]
    else:
        return slice_


# TODO: random_axes="v", axis_order="vtcz" processes whole FOVs time-series in random FOV order
# (better for testing)
def send_nd2(filename, axis_order="tvcz", slices={}, delayed=True):
    if delayed is True:
        delayed = dask.delayed(pure=True)
    elif delayed is False:
        delayed = lambda func, **kwargs: func
    nd2 = nd2reader.ND2Reader(filename)
    iterators = [
        _select_indices(
            np.arange(nd2.sizes.get(axis_name, 1)), slices.get(axis_name, slice(None))
        )
        for axis_name in axis_order
    ]
    iterator = product(*iterators)
    # send whole-file metadata
    # TODO: we don't wrap in delayed, should we?
    nd2_metadata = parse_nd2_metadata(nd2)
    yield {"type": "nd2_metadata", "metadata": nd2_metadata}
    # send frames
    for idxs in iterator:
        coords = dict(zip(axis_order, idxs))
        # TODO: this might not work delayed because ND2Reader doesn't reopen file handles correctly
        # TODO: use an LRU buffer for open ND2 file handles?
        # image = delayed(nd2.get_frame_2D)(**coords)
        fov_num = coords["v"]
        channel = nd2.metadata["channels"][coords["c"]]
        z_level = coords["z"]
        t = coords["t"]
        image = delayed(get_nd2_frame)(filename, fov_num, coords["c"], t)
        # we can put arbitrary per-frame metadata here
        # TODO: do we need a way to encode metadata that differs between individual frames? [maybe not.]
        image_metadata = {
            "dummy_metadata": 0,
            "channel": channel,
            "z_level": z_level,
            "fov_num": fov_num,
            "t": t,
        }
        # TODO: we will replace the above with a more elegant way of encoding
        # channel, timepoint, FOV number, z-level, etc.
        msg = {"type": "image", "image": image, "metadata": image_metadata}
        yield msg


def _slice_directory_keys(keys, slices):
    idxs = tuple(map(np.unique, zip(*keys)))
    selected_idxs = tuple(map(_select_indices, idxs, slices))
    for key in sorted(keys):
        for idx, selected_idx in zip(key, selected_idxs):
            if idx not in selected_idx:
                break
        else:
            yield key
    return selected_idxs


def slice_directory(root_dir, pattern, axis_order="tvc", slices={}):
    pattern = re.compile(pattern)
    root_dir = Path(root_dir)
    unmatched_filenames = []
    keys_to_filenames = {}
    for file in root_dir.iterdir():
        match = pattern.match(file.name)
        if match:
            idxs = {
                k: excepts(ValueError, int, lambda _: v)(v)
                for k, v in match.groupdict().items()
            }
            key = tuple([idxs[axis] for axis in axis_order])
            keys_to_filenames[key] = file.name
        else:
            unmatched_filenames.append(file.name)
    selected_keys = _slice_directory_keys(
        keys_to_filenames.keys(), [slices.get(axis, slice(None)) for axis in axis_order]
    )
    iterator = ((key, root_dir / keys_to_filenames[key]) for key in selected_keys)
    return iterator, unmatched_filenames


def get_eaton_fish_frame(filename):
    with h5py.File(filename) as f:
        # load array into memory as a numpy array
        # TODO: do we ever want to memmap?
        frame = f["data"][()]
    return frame


def send_eaton_fish(root_dir, pattern, axis_order="tvc", slices={}, delayed=True):
    root_dir = Path(root_dir)
    if delayed is True:
        delayed = dask.delayed(pure=True)
    elif delayed is False:
        delayed = lambda func, **kwargs: func
    filenames, unmatched_filenames = slice_directory(
        root_dir, pattern, axis_order=axis_order, slices=slices
    )
    for suffix in ("initial.hdf5", "init_fixation.hdf5", "fixed.hdf5"):
        filename = root_dir / suffix
        if not filename.exists():
            continue
        image = delayed(get_eaton_fish_frame)(filename)
        image_metadata = {}  # TODO
        msg = {
            "type": "image",
            "image": image,
            "metadata": image_metadata,
            "image_type": "eaton_fish_fixation",
            "image_name": filename.stem,
        }
        yield msg
    for key, filename in filenames:
        image = delayed(get_eaton_fish_frame)(filename)
        coords = dict(zip(axis_order, key))
        image_metadata = {
            "dummy_metadata": 0,
            "channel": coords["c"],
            "fov_num": coords["v"],
            "t": coords["t"],
        }
        msg = {
            "type": "image",
            "image": image,
            "metadata": image_metadata,
            "image_type": "eaton_fish",
        }
        yield msg
    # TODO: these are not sent sorted, but it shouldn't matter
    for metadata_filename in root_dir.glob("metadata_*.hdf5"):
        metadata = delayed(pd.read_hdf)(metadata_filename)
        msg = {"type": "table", "table": metadata, "table_type": "eaton_fish_metadata"}
        yield msg


def get_hdf5_frame(filename, hdf5_path, slice_=()):
    return h5py.File(filename)[hdf5_path][slice_]


def send_hdf5(filename, delayed=True):
    if delayed is True:
        delayed = dask.delayed(pure=True)
    elif delayed is False:
        delayed = lambda func, **kwargs: func
    h5 = h5py.File(filename)
    # TODO: slicing
    # TODO: send whole-file metadata
    for channel, channel_group in h5.items():
        # TODO: fov -> fov_num?, z -> z_level (also in convert_nd2_to_hdf5)
        for fov, ary in channel_group.items():
            for z in range(ary.shape[0]):
                for t in range(ary.shape[1]):
                    image = delayed(get_hdf5_frame)(
                        filename, f"{channel}/{fov}", (z, t)
                    )
                    image_metadata = {
                        "dummy_metadata": 0,
                        "channel": channel,
                        "z_level": z,
                        "fov_num": fov,
                        "t": t,
                    }
                    msg = {"type": "image", "image": image, "metadata": image_metadata}
                    yield msg


def convert_nd2_to_hdf5(
    nd2_filename,
    hdf5_path,
    file_axes=["t", "fov"],
    dataset_axes=["channel"],
    slice_axes=None,
    chunks=dict(
        z=None,
        t=5,
        channel_idx=1,
        fov=1,
        height=None,
        width=None,
    ),
    slices={},
):
    file_and_dataset_axes = set(file_axes) | set(dataset_axes)
    unknown_axes = file_and_dataset_axes - set(
        ["t", "fov", "channel", "channel_num", "z"]
    )
    if unknown_axes:
        raise ValueError(f"unknown axes: {', '.join(unknown_axes)}")
    if slice_axes is None:
        slice_axes = []
        if not set(["fov", "fov_idx"]) & file_and_dataset_axes:
            slice_axes.append("fov_idx")
        if not set(["t", "t_idx"]) & file_and_dataset_axes:
            slice_axes.append("t_idx")
        if not set(["channel", "channel_num"]) & file_and_dataset_axes:
            slice_axes.append("channel_idx")
        if "z" not in file_and_dataset_axes:
            slice_axes.append("z_idx")
    # the keys given to dataset_shape_func and dataset_chunks_func
    # have keys named channel, t (not channel_idx, t_idx)
    dataset_creation_axes = [axis.split("_")[0] for axis in slice_axes]
    hdf5_filename = str(hdf5_path)
    if file_axes:
        hdf5_filename += "_".join([f"{a}={{{a}}}" for a in file_axes])
    if not (
        hdf5_filename.lower().endswith(".hdf5") or hdf5_filename.lower().endswith(".h5")
    ):
        hdf5_filename += r".hdf5"
    filename_func = lambda key: hdf5_filename.format(
        **{axis: v for axis, v in key.items() if axis in file_axes}
    )
    dataset_func = lambda key: "/".join(f"{axis}={key[axis]}" for axis in dataset_axes)
    slice_func = lambda key: (
        *(key[axis] for axis in slice_axes),
        slice(None),
        slice(None),
    )
    dataset_shape_func = lambda key: (
        *(len(key[axis]) for axis in dataset_creation_axes),
        key["height"],
        key["width"],
    )
    dataset_chunks_func = lambda key: tuple(
        dataset_shape_func(key)[axis_idx]
        if chunks[axis] is None
        else min(chunks[axis], dataset_shape_func(key)[axis_idx])
        for axis_idx, axis in enumerate([*dataset_creation_axes, "height", "width"])
    )
    _convert_nd2_to_hdf5(
        nd2_filename,
        filename_func,
        dataset_func,
        slice_func,
        dataset_shape_func,
        dataset_chunks_func,
        slices=slices,
    )


def _convert_nd2_to_hdf5(
    nd2_filename,
    filename_func,
    dataset_func,
    slice_func,
    dataset_shape_func,
    dataset_chunks_func,
    slices={},
):
    if "channel" in slices and "channel_num" in slices:
        raise ValueError("cannot specify both channel and channel_num slices")
    nd2 = nd2reader.ND2Reader(nd2_filename)
    nd2_channels = nd2.metadata["channels"]
    frame = nd2.get_frame_2D()
    dtype = frame.dtype
    shape = frame.shape
    del frame
    if "channel" in slices:
        channel_slices = slices["channel"]
        if isinstance(channel_slices, slice):
            slices["channel_num"] = slice(
                nd2_channels.index(channel_slices.start),
                nd2_channels.index(channel_slices.stop),
                channel_slices.step,
            )
        else:
            slices["channel_num"] = [nd2_channels.index(c) for c in channel_slices]
    channel_nums = _select_indices(
        np.arange(nd2.sizes.get("c", 1)), slices.get("channel_num", slice(None))
    )
    fovs = _select_indices(
        np.arange(nd2.sizes.get("v", 1)), slices.get("fov", slice(None))
    )
    zs = _select_indices(np.arange(nd2.sizes.get("z", 1)), slices.get("z", slice(None)))
    ts = _select_indices(np.arange(nd2.sizes.get("t", 1)), slices.get("t", slice(None)))
    hdf5_files = {}
    try:
        for channel_idx, channel_num in enumerate(
            tqdm(channel_nums, desc="c", leave=None)
        ):
            channel = nd2_channels[channel_idx]
            for fov_idx, fov in enumerate(tqdm(fovs, desc="v", leave=None)):
                for z_idx, z in enumerate(tqdm(zs, desc="z", leave=None)):
                    for t_idx, t in enumerate(tqdm(ts, desc="t", leave=None)):
                        key = dict(
                            fov=fov,
                            fov_idx=fov_idx,
                            channel=channel,
                            channel_num=channel_num,
                            channel_idx=channel_idx,
                            z=z,
                            z_idx=z_idx,
                            t=t,
                            t_idx=t_idx,
                        )
                        frame = nd2.get_frame_2D(c=channel_num, v=fov, z=z, t=t)
                        hdf5_filename = filename_func(key)
                        hdf5_file = hdf5_files.get(hdf5_filename)
                        if hdf5_file is None:
                            hdf5_file = hdf5_files[hdf5_filename] = h5py.File(
                                hdf5_filename, "a"
                            )
                        dataset_path = dataset_func(key)
                        if dataset_path not in hdf5_file:
                            dataset_key = dict(
                                # this key is called channel, not channel_nums
                                # for naming consistency in the key that is
                                # passed to dataset_shape_func and
                                # dataset_chunks_func
                                channel=channel_nums,
                                fov=fovs,
                                z=zs,
                                t=ts,
                                height=shape[0],
                                width=shape[1],
                            )
                            dataset_shape = dataset_shape_func(dataset_key)
                            dataset_chunks = dataset_chunks_func(dataset_key)
                            hdf5_file.create_dataset(
                                dataset_path,
                                shape=dataset_shape,
                                chunks=dataset_chunks,
                                dtype=dtype,
                            )
                        slice_ = slice_func(key)
                        hdf5_file[dataset_path][slice_] = frame

    finally:
        for hdf5_file in hdf5_files.values():
            hdf5_file.close()
