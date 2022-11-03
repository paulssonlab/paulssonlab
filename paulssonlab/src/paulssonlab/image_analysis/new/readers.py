import numpy as np
import dask
import h5py
import nd2reader
import cachetools
from tqdm.auto import tqdm
from cytoolz import excepts
from itertools import product
from paulssonlab.io.metadata import parse_nd2_metadata

ND2READER_CACHE = cachetools.LFUCache(maxsize=48)


def _get_nd2_reader(filename, **kwargs):
    # TODO: removed memmap flag
    return nd2reader.ND2Reader(filename, **kwargs)


get_nd2_reader = cachetools.cached(cache=ND2READER_CACHE)(_get_nd2_reader)


def get_nd2_frame(filename, position, channel_idx, t):
    reader = get_nd2_reader(filename)
    ary = reader.get_frame_2D(v=position, c=channel_idx, t=t)
    return ary


def _select_indices(idxs, slice_):
    if isinstance(slice_, slice):
        return idxs[slice_]
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
        frame = delayed(get_eaton_fish_frame)(filename)
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


def convert_nd2_to_hdf5(nd2_filename, hdf5_filename, slices):
    nd2 = nd2reader.ND2Reader(nd2_filename)
    frame = nd2.get_frame_2D()
    dtype = frame.dtype
    shape = frame.shape
    del frame
    if "c" in slices and not isinstance(slices["c"], slice):
        slices["c"] = [nd2.metadata["channels"].index(c) for c in slices["c"]]
    channel_idxs = _select_indices(
        np.arange(nd2.sizes.get("c", 1)), slices.get("c", slice(None))
    )
    fovs = _select_indices(
        np.arange(nd2.sizes.get("v", 1)), slices.get("v", slice(None))
    )
    zs = _select_indices(np.arange(nd2.sizes.get("z", 1)), slices.get("z", slice(None)))
    ts = _select_indices(np.arange(nd2.sizes.get("t", 1)), slices.get("t", slice(None)))
    with h5py.File(hdf5_filename, "a") as f:
        for channel_idx in tqdm(channel_idxs, desc="c", leave=None):
            channel = nd2.metadata["channels"][channel_idx]
            for fov in tqdm(fovs, desc="v", leave=None):
                ary = f.create_dataset(
                    f"{channel}/{fov}",
                    shape=(len(zs), len(ts), *shape),
                    chunks=(1, 5, *shape),
                    dtype=dtype,
                )
                for z_idx, z in enumerate(tqdm(zs, desc="z", leave=None)):
                    for t_idx, t in enumerate(tqdm(ts, desc="t", leave=None)):
                        ary[z_idx, t_idx, :, :] = nd2.get_frame_2D(
                            c=channel_idx, v=fov, z=z, t=t
                        )
