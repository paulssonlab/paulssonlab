import numpy as np
import nd2reader
from numcodecs import Blosc, Delta
import operator
from functools import partial, reduce
import itertools
from more_itertools import rstrip
from collections.abc import Sequence, Iterable
from copy import deepcopy
from utils import tqdm_auto, open_zarr_group

DEFAULT_FRAME_COMPRESSOR = Blosc(
    cname="zstd", clevel=5, shuffle=Blosc.SHUFFLE, blocksize=0
)
DEFAULT_FRAME_CHUNKS = (1, 1, 512, 512)


def map_ndarray(func, ary, axis=0, out=None, progress_bar=tqdm_auto):
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
    progress_bar=tqdm_auto,
):
    res = None
    shape = None
    dtype = None
    for pos_name, position in progress_bar(in_group.items()):
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
        if res is None:
            map_ndarray(func, position, axis=axis, out=ary, progress_bar=progress_bar)
        else:
            ary[:] = res
            res = None
        ary.attrs["processed"] = True
    return out_group


def process_attrs(func, in_group, out_group, progress_bar=tqdm_auto):
    for pos_name, position in progress_bar(in_group.items()):
        pos_group = out_group.require_group(pos_name)
        if "processed" and pos_group.attrs and pos_group.attrs["processed"]:
            continue
        else:
            pos_group.attrs["processed"] = False
        res = func(position)
        pos_group["result"] = res
        pos_group.attrs["processed"] = True
    return out_group


def ingest_nd2_file(nd2_path, zarr_path, **kwargs):
    nd2 = nd2reader.ND2Reader(nd2_path)
    raw_group = open_zarr_group(zarr_path).require_group("raw")
    return ingest_nd2(nd2, raw_group, **kwargs)


def ingest_nd2(
    nd2,
    raw_group,
    progress_bar=tqdm_auto,
    compressor=DEFAULT_FRAME_COMPRESSOR,
    chunks=DEFAULT_FRAME_CHUNKS,
    **kwargs,
):
    meta = deepcopy(nd2.metadata)
    meta["date"] = meta["date"].isoformat()
    raw_group.attrs["metadata"] = meta
    for v in tqdm_auto(range(nd2.sizes["v"])):
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
        ary.attrs["metadata"] = meta
        for c in progress_bar(range(nd2.sizes["c"]), leave=False):
            for t in progress_bar(range(nd2.sizes["t"]), leave=False):
                ary[c, t, :, :] = nd2.get_frame_2D(c=c, t=t, v=v)
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


def quantize_frames(in_group, out_group, bits, random=True, progress_bar=tqdm_auto):
    out_group.attrs["quantization"] = {"bits": bits, "random": random}
    return process_positions(
        lambda frame: quantize_frame(frame, bits, random=random),
        in_group,
        out_group,
        axis=(0),
    )
