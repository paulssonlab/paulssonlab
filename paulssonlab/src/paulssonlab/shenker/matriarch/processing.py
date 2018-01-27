import numpy as np
import nd2reader
from numcodecs import Blosc, Delta
import operator
from functools import partial, reduce
import itertools
from more_itertools import rstrip
from collections.abc import Sequence, Iterable
from copy import deepcopy
import time
from natsort import natsorted
from utils import tqdm_auto, open_zarr_group, timestamp_to_isoformat

DEFAULT_FRAME_COMPRESSOR = Blosc(
    cname="zstd", clevel=5, shuffle=Blosc.SHUFFLE, blocksize=0
)
DEFAULT_FRAME_CHUNKS = (1, 1, 512, 512)
DEFAULT_PROGRESS_BAR = partial(tqdm_auto, smoothing=0.8)


def map_ndarray(func, ary, axis=0, out=None, progress_bar=DEFAULT_PROGRESS_BAR):
    if not isinstance(axis, Sequence):
        axis = [axis]
    iter_shape = [ary.shape[dim] for dim in axis]
    _iter = itertools.product(*[range(l) for l in iter_shape])
    if progress_bar is not None:
        _iter = progress_bar(_iter, total=reduce(operator.mul, iter_shape))
    for iter_idxs in _iter:
        idxs = [slice(None)] * len(ary.shape)
        for dim, iter_idx in zip(axis, iter_idxs):
            idxs[dim] = iter_idx
        idxs = tuple(rstrip(idxs, lambda x: x == slice(None)))
        if out is None:
            new_ary = func(ary[idxs])
            iter_shape = tuple(map(lambda dim: ary.shape[dim], axis))
            new_shape = iter_shape + new_ary.shape
            out = np.zeros(new_shape, dtype=new_ary.dtype)
            out[idxs] = new_ary
        else:
            out[idxs] = func(ary[idxs])
    return out


def process_positions(
    func,
    in_group,
    out_group,
    axis=0,
    compressor=None,
    chunks=None,
    order=None,
    progress_bar=DEFAULT_PROGRESS_BAR,
):
    res = None
    shape = None
    dtype = None
    pbar = progress_bar(natsorted(in_group.items(), key=operator.itemgetter(0)))
    for pos_name, position in pbar:
        pbar.set_description("position {}".format(pos_name))
        time_started = time.time()
        if compressor is None:
            compressor = position.compressor
        if order is None:
            order = position.order
        if pos_name in out_group:
            ary = out_group[pos_name]
            shape = ary.shape
            dtype = ary.dtype
            if "processed" and ary.attrs and ary.attrs["processed"]:
                continue
        else:
            if shape is None or dtype is None:
                res = map_ndarray(
                    func,
                    position,
                    axis=axis,
                    progress_bar=partial(progress_bar, leave=False),
                )
                shape = res.shape
                dtype = res.dtype
            if chunks is None and len(shape) == len(position.shape):
                chunks = position.chunks
            ary = out_group.zeros(
                pos_name,
                shape=shape,
                dtype=dtype,
                compressor=compressor,
                chunks=chunks,
                order=order,
            )
        ary.attrs["processed"] = False
        ary.attrs["time_started"] = timestamp_to_isoformat(time_started)
        if res is None:
            map_ndarray(func, position, axis=axis, out=ary, progress_bar=progress_bar)
        else:
            ary[:] = res
            res = None
        time_finished = time.time()
        ary.attrs["time_finished"] = timestamp_to_isoformat(time_finished)
        ary.attrs["time_elapsed"] = time_finished - time_started
        ary.attrs["processed"] = True
    return out_group


def process_attrs(func, in_group, out_group, progress_bar=DEFAULT_PROGRESS_BAR):
    for pos_name, position in progress_bar(in_group.items()):
        pos_group = out_group.require_group(pos_name)
        if "processed" and pos_group.attrs and pos_group.attrs["processed"]:
            continue
        pos_group.attrs["processed"] = False
        time_started = time.time()
        pos_group.attrs["time_started"] = timestamp_to_isoformat(time_started)
        res = func(position)
        pos_group.attrs["result"] = res
        time_finished = time.time()
        pos_group.attrs["time_finished"] = timestamp_to_isoformat(time_finished)
        pos_group.attrs["time_elapsed"] = time_finished - time_started
        pos_group.attrs["processed"] = True
    return out_group


def ingest_nd2_file(nd2_path, zarr_path, **kwargs):
    nd2 = nd2reader.ND2Reader(nd2_path)
    raw_group = open_zarr_group(zarr_path).require_group("raw")
    return ingest_nd2(nd2, raw_group, **kwargs)


def ingest_nd2(
    nd2,
    raw_group,
    progress_bar=DEFAULT_PROGRESS_BAR,
    compressor=DEFAULT_FRAME_COMPRESSOR,
    chunks=DEFAULT_FRAME_CHUNKS,
    **kwargs,
):
    meta = deepcopy(nd2.metadata)
    meta["date"] = meta["date"].isoformat()
    raw_group.attrs["metadata"] = meta
    pbar_v = progress_bar(range(nd2.sizes["v"]))
    for v in pbar_v:
        pbar_v.set_description("position {}".format(v))
        ary = raw_group.require_dataset(
            "{:d}".format(v),
            shape=[nd2.sizes[n] for n in "ctyx"],
            chunks=chunks,
            dtype="u2",
            order="C",
            compressor=compressor,
            **kwargs,
        )
        if "ingested" in ary.attrs and ary.attrs["ingested"]:
            continue
        ary.attrs["ingested"] = False
        time_started = time.time()
        ary.attrs["time_started"] = timestamp_to_isoformat(time_started)
        ary.attrs["metadata"] = meta
        pbar_c = progress_bar(range(nd2.sizes["c"]), leave=False)
        for c in pbar_c:
            pbar_c.set_description("channel {}".format(c))
            pbar_t = progress_bar(range(nd2.sizes["t"]), leave=False)
            for t in pbar_t:
                pbar_t.set_description("timepoint {}".format(t))
                ary[c, t, :, :] = nd2.get_frame_2D(c=c, t=t, v=v)
        time_finished = time.time()
        ary.attrs["time_finished"] = timestamp_to_isoformat(time_finished)
        ary.attrs["time_elapsed"] = time_finished - time_started
        ary.attrs["ingested"] = True
    return raw_group


def quantize_frame(arr, bits, random=True):
    factor = 2**bits
    if random:
        return np.floor((arr / factor) + np.random.random(size=arr.shape)).astype(
            arr.dtype
        )
    else:
        return arr // factor


def quantize_frames(
    in_group, out_group, bits, random=True, progress_bar=DEFAULT_PROGRESS_BAR
):
    out_group.attrs["quantization"] = {"bits": bits, "random": random}
    out_group.attrs["metadata"] = in_group.attrs["metadata"]
    return process_positions(
        lambda frame: quantize_frame(frame, bits, random=random),
        in_group,
        out_group,
        axis=(0),
    )
