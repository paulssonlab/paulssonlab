import numpy as np
import pandas as pd
import holoviews as hv
import dask
from dask import delayed
import dask.array as da
import nd2reader
import skimage.filters
import skimage.feature
from skimage.feature import hessian_matrix, hessian_matrix_eigvals
import scipy.ndimage
from cytoolz import reduce, compose
import numpy_indexed
import cachetools
from numbers import Integral
import warnings

ND2READER_CACHE = cachetools.LFUCache(maxsize=48)


def _get_nd2_reader(filename, memmap=False, **kwargs):
    return nd2reader.ND2Reader(filename, memmap=memmap, **kwargs)


get_nd2_reader = cachetools.cached(cache=ND2READER_CACHE)(_get_nd2_reader)
# get_nd2_reader = _get_nd2_reader

# TODO: modified
def get_nd2_frame(filename, position, channel, t, memmap=False):
    reader = get_nd2_reader(filename, memmap=memmap)
    if not isinstance(channel, Integral):
        channel = reader.metadata["channels"].index(channel)
    ary = reader.get_frame_2D(v=position, c=channel, t=t, memmap=memmap)
    return ary


# TODO: new
def nd2_to_dask(filename, position, channel):
    nd2 = get_nd2_reader(filename)
    frame0 = get_nd2_frame(filename, position, channel, 0)
    frames = [
        delayed(get_nd2_frame)(filename, position, channel, t)
        for t in range(nd2.sizes["t"])
    ]
    arrays = [
        da.from_delayed(frame, dtype=frame0.dtype, shape=frame0.shape)
        for frame in frames
    ]
    stack = da.stack(arrays, axis=0)
    return stack


# TODO: new
def nd2_to_futures(client, filename, position, channel):
    nd2 = get_nd2_reader(filename)
    frames = [
        client.submit(get_nd2_frame, filename, position, channel, t)
        for t in range(nd2.sizes["t"])
    ]
    return frames


def RevImage(img, **kwargs):
    return _RevImage(hv.Image, img, **kwargs)


def _RevImage(cls, img, **kwargs):
    return cls(img[::-1], bounds=(0, 0, img.shape[1], img.shape[0])).options(
        invert_yaxis=True
    )


# TODO: changed argument order, use in matriarch measure func
def map_over_labels(func, label_image, intensity_image):
    # assumes are consecutive integers 0,...,N
    groups = numpy_indexed.group_by(
        label_image.ravel(), intensity_image.ravel(), reduction=func
    )
    return [g[1] for g in groups]


def repeat_apply(func, n):
    if n <= 0:
        return lambda x: x
    return reduce(lambda f1, f2: compose(f1, f2), [func] * n)


def hessian_eigenvalues(img):
    return hessian_matrix_eigvals(hessian_matrix(img.astype(np.float32), order="rc"))


def gaussian_box_approximation(ary, sigma, n=3, mode="nearest"):
    w_ideal = np.sqrt((12 * sigma**2 / n) + 1)
    w_l = int(np.floor(w_ideal))
    if w_l % 2 == 0:
        w_l -= 1
    w_u = w_l + 2
    m = (12 * sigma**2 - n * w_l**2 - 4 * n * w_l - 3 * n) / (-4 * w_l - 4)
    m = round(m)
    for i in range(m):
        ary = scipy.ndimage.filters.uniform_filter(ary, w_l, mode=mode)
    for i in range(n - m):
        ary = scipy.ndimage.filters.uniform_filter(ary, w_u, mode=mode)
    return ary


def normalize_componentwise(
    img,
    img_labels,
    weighted=False,
    label_index=None,
    dilation=5,
    in_place=False,
    dtype=np.float32,
):
    if not in_place:
        img = img.astype(dtype).copy()
    if label_index is None:
        label_index = np.unique(img_labels)
    maxes = scipy.ndimage.maximum(img, labels=img_labels, index=label_index)
    if weighted:
        median = np.median(maxes[1:])  # exclude background from median
    img_labels = repeat_apply(skimage.morphology.dilation, dilation)(img_labels)
    img[img_labels == 0] = 0
    for idx, label in enumerate(label_index):
        mask = img_labels == label
        if weighted:
            img[mask] *= min(maxes[idx] / median, 1)
        else:
            img[mask] /= maxes[idx]
    return img


DEFAULT_REGIONPROPS = [
    "area",
    "centroid",
    "eccentricity",
    "min_intensity",
    "mean_intensity",
    "max_intensity",
    "major_axis_length",
    "minor_axis_length",
    "orientation",
    "perimeter",
    "solidity",
    "weighted_centroid",
]


def get_regionprops(label_image, intensity_image, properties=DEFAULT_REGIONPROPS):
    rps = skimage.measure.regionprops(
        label_image, intensity_image, coordinates="rc", cache=False
    )
    if not len(rps):
        return None
    cols = {prop: [getattr(rp, prop) for rp in rps] for prop in properties}
    for col, values in list(cols.items()):
        if isinstance(values[0], tuple):
            del cols[col]
            # TODO: store coordinates as multiindex?
            cols[col + "_x"] = [v[0] for v in values]
            cols[col + "_y"] = [v[1] for v in values]
    df = pd.DataFrame(cols, index=range(1, len(rps) + 1))
    df.index.name = "label"
    return df
