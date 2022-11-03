import numpy as np
import dask
import h5py
import nd2reader
import cachetools
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


# TODO: random_axes="v", axis_order="vtcz" processes whole FOVs time-series in random FOV order
# (better for testing)
def send_nd2(filename, axis_order="tvcz", slices={}, delayed=True):
    if delayed is True:
        delayed = dask.delayed(pure=True)
    elif delayed is False:
        delayed = lambda func, **kwargs: func
    nd2 = nd2reader.ND2Reader(filename)
    iterators = [
        np.arange(nd2.sizes.get(axis_name, 1))[slices.get(axis_name, slice(None))]
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
    def _select_indices(idxs, slice_):
        if isinstance(slice_, slice):
            return idxs[slice_]
        else:
            return slice_

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
