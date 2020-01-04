import numpy as np
import pandas as pd
import holoviews as hv
import dask
from dask import delayed
import dask.array as da
import nd2reader
import skimage
import skimage.filters
import skimage.feature
from skimage.feature import hessian_matrix, hessian_matrix_eigvals
import scipy.ndimage
from cytoolz import reduce, compose, partial
import numpy_indexed
import cachetools
from numbers import Integral
import warnings
import zarr
from numcodecs import Blosc
from collections import Sequence, Mapping, defaultdict


def tree():
    return defaultdict(tree)


DEFAULT_COMPRESSOR = Blosc(cname="zstd", clevel=5, shuffle=Blosc.SHUFFLE, blocksize=0)
DEFAULT_ORDER = "C"
zarrify = partial(
    zarr.array, compressor=DEFAULT_COMPRESSOR, order=DEFAULT_ORDER, chunks=False
)

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


def RevImage(img, **kwargs):
    return _RevImage(hv.Image, img, **kwargs)


def _RevImage(cls, img, **kwargs):
    return cls(img[::-1], bounds=(0, 0, img.shape[1], img.shape[0])).options(
        invert_yaxis=True
    )


# TODO: new. incorporate into recursive_map?
def recursive_sequence_map(func, data, max_level=0):
    if max_level == 0:
        return func(data)
    else:
        values = type(data)(
            [recursive_sequence_map(func, d, max_level=max_level - 1) for d in data]
        )
        values = func(values)
        return values


# FROM: https://stackoverflow.com/questions/42095393/python-map-a-function-over-recursive-iterables/42095505
# TODO: document!!! and replace map_collections
# TODO: incorporate recursive_sequence_map?
def recursive_map(
    func,
    data,
    shortcircuit=(),
    ignore=(),
    keys=False,
    max_level=None,
    predicate=None,
    key_predicate=None,
):
    if max_level is not None and max_level is not False:
        if max_level == 0:
            return func(data)
        max_level -= 1
    apply = lambda x: recursive_map(
        func,
        x,
        shortcircuit=shortcircuit,
        ignore=ignore,
        keys=keys,
        max_level=max_level,
        predicate=predicate,
        key_predicate=key_predicate,
    )
    if isinstance(data, shortcircuit):
        return func(data)
    # to avoid an infinite recursion error
    # we need to hardcode that strs are ignored, because str[0] is a str, and hence a Sequence
    elif isinstance(data, str) or ignore is not True and isinstance(data, ignore):
        return data
    elif isinstance(data, Mapping):
        if keys:
            values = {apply(k): apply(v) for k, v in data.items()}
        else:
            values = {k: apply(v) for k, v in data.items()}
        if key_predicate is not None:
            values = {k: v for k, v in values.items() if key_predicate(k)}
        if predicate is not None:
            values = {k: v for k, v in values.items() if predicate(v)}
        return type(data)(values)
    elif isinstance(data, Sequence):
        values = [apply(v) for v in data]
        if predicate is not None:
            values = [v for v in values if predicate(v)]
        return type(data)(values)
    elif ignore is True or True in ignore:
        return data
    else:
        return func(data)


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


def hessian_eigenvalues(img, sigma=1.5):
    return hessian_matrix_eigvals(
        hessian_matrix(skimage.img_as_float(img), sigma, mode="nearest", order="rc")
    )


# def hessian_eigenvalues(img, sigma=1.5):
#     I = skimage.filters.gaussian(img, sigma)
#     I_x = skimage.filters.sobel_h(I)
#     I_y = skimage.filters.sobel_v(I)
#     I_xx = skimage.filters.sobel_h(I_x)
#     I_xy = skimage.filters.sobel_v(I_x)
#     I_yx = skimage.filters.sobel_h(I_y)
#     I_yy = skimage.filters.sobel_v(I_y)
#     kappa_1 = (I_xx + I_yy) / 2
#     with warnings.catch_warnings():
#         warnings.simplefilter('ignore', RuntimeWarning)
#         kappa_2 = (np.sqrt((I_xx + I_yy)**2 - 4*(I_xx*I_yy - I_xy*I_yx))) / 2
#     k1 = kappa_1 + kappa_2
#     k2 = kappa_1 - kappa_2
#     k1[np.isnan(k1)] = 0
#     k2[np.isnan(k2)] = 0
#     return k1, k2

# TODO: modified
def gaussian_box_approximation(ary, sigma, n=3, mode="nearest", cval=0):
    w_ideal = np.sqrt((12 * sigma**2 / n) + 1)
    w_l = int(np.floor(w_ideal))
    if w_l % 2 == 0:
        w_l -= 1
    w_u = w_l + 2
    m = (12 * sigma**2 - n * w_l**2 - 4 * n * w_l - 3 * n) / (-4 * w_l - 4)
    m = round(m)
    for i in range(m):
        ary = scipy.ndimage.filters.uniform_filter(ary, w_l, mode=mode, cval=0)
    for i in range(n - m):
        ary = scipy.ndimage.filters.uniform_filter(ary, w_u, mode=mode, cval=0)
    return ary


# TODO: modified
def normalize_componentwise(
    img, img_labels, label_index=None, dilation=10, in_place=False
):
    img = skimage.img_as_float(img, force_copy=(not in_place))
    if label_index is None:
        label_index = np.unique(img_labels)
    img_labels = repeat_apply(skimage.morphology.dilation, dilation)(img_labels)
    mins, maxes, _, _ = scipy.ndimage.extrema(img, labels=img_labels, index=label_index)
    img[img_labels == 0] = 0
    for idx, label in enumerate(label_index):
        if label == 0:
            continue
        mask = img_labels == label
        img[mask] = (img[mask] - mins[idx]) / (maxes[idx] - mins[idx])
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
