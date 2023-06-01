import warnings

import numba
import numpy as np
import pandas as pd
import scipy
import skimage.morphology
from cytoolz import compose
from skimage.feature import hessian_matrix, hessian_matrix_eigvals

from paulssonlab.image_analysis.blur import scipy_box_blur
from paulssonlab.image_analysis.util import repeat_apply

gaussian_box_blur = scipy_box_blur
del scipy_box_blur


def quantize(img, bits, random=True):
    factor = 2**bits
    if random:
        return np.floor((img / factor) + np.random.random(size=img.shape)).astype(
            img.dtype
        )
    else:
        return img // factor


def downsample(img, factor):
    return img[::factor, ::factor]


def _accumulator_dtype(dtype):
    if np.issubdtype(dtype, np.bool_):
        return np.uint64
    elif np.issubdtype(dtype, np.signedinteger):
        return np.int64
    elif np.issubdtype(dtype, np.unsignedinteger):
        return np.uint64
    elif np.issubdtype(dtype, np.floating):
        return np.float64
    elif np.issubdtype(dtype, np.complex):
        return np.complex128
    else:
        return NotImplementedError


# FROM: https://alyssaq.github.io/2014/understanding-hough-transform/
def hough_line_intensity(img, theta=None):
    if theta is None:
        theta = np.linspace(-np.pi / 2, np.pi / 2, 180)
    width, height = img.shape
    diagonal = int(np.ceil(np.sqrt(width**2 + height**2)))
    rho = np.linspace(-diagonal, diagonal, diagonal * 2)
    accumulator = np.zeros(
        (2 * diagonal, len(theta)), dtype=_accumulator_dtype(img.dtype)
    )
    _hough_line_intensity(accumulator, img, theta, diagonal)
    return accumulator, theta, rho


# TODO: fastmath doesn't seem to do anything; could also try parallel=True
@numba.njit
def _hough_line_intensity(accumulator, img, theta, diagonal):
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    num_thetas = len(theta)
    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            if img[i, j] == 0:  # optimization for boolean input images
                continue
            for k in range(num_thetas):
                # TODO: is round correct here?
                rho = int(np.round_(i * sin_theta[k] + j * cos_theta[k])) + diagonal
                accumulator[rho, k] += img[i, j]


# TODO: take out kwarg instead of in_place, to conform with skimage.morphology.remove_small_objects
def remove_large_objects(labeled_img, max_size, in_place=False):
    if not in_place:
        labeled_img = labeled_img.copy()
    component_sizes = np.bincount(labeled_img.ravel())
    too_big = component_sizes > max_size
    too_big_mask = too_big[labeled_img]
    labeled_img[too_big_mask] = 0
    return labeled_img


# TODO: I think weighted=True is useless, it might do the opposite of what you want...
def normalize_componentwise(
    img,
    img_labels,
    ensure_nonnegative=True,
    weighted=False,
    label_index=None,
    inverse_label_map=None,
    dilation=5,
    in_place=False,
    dtype=np.float32,
):
    if not in_place:
        img = img.astype(dtype).copy()
    if ensure_nonnegative:
        img -= np.nanmin(img)
    if label_index is None:
        label_index = np.unique(img_labels)
        # remove background
        if label_index[0] == 0:
            label_index = label_index[1:]
    if inverse_label_map is None:
        inverse_label_map = np.zeros(label_index.max(), dtype=np.uint16)
        for idx, label in enumerate(label_index):
            inverse_label_map[label - 1] = idx
    maxes = scipy.ndimage.maximum(img, labels=img_labels, index=label_index)
    if weighted:
        median = np.nanmedian(maxes[~np.isinf(maxes)])
        # if any scales are inf, that means that there are no pixels with that label
        # so won't matter in _normalize_componentwise
        scale = np.minimum(maxes / median, 1)
    else:
        with np.errstate(divide="ignore"):
            scale = 1 / maxes
    img_labels = repeat_apply(skimage.morphology.dilation, dilation)(img_labels)
    img[img_labels == 0] = 0
    return _normalize_componentwise(img, img_labels, inverse_label_map, scale)


@numba.njit
def _normalize_componentwise(
    img,
    img_labels,
    inverse_label_map,
    scale,
):
    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            if img_labels[i, j] == 0:
                continue
            img[i, j] *= scale[inverse_label_map[img_labels[i, j] - 1]]
    return img


# FROM: https://stackoverflow.com/questions/21242011/most-efficient-way-to-calculate-radial-profile
def radial_profile(img, center=None):
    if center is None:
        center = (img.shape[0] // 2, img.shape[1] // 2)
    y, x = np.indices((img.shape))
    r = np.sqrt((x - center[0]) ** 2 + (y - center[1]) ** 2)
    r = r.astype(np.int_)
    counts = np.bincount(r.ravel(), img.ravel())
    bin_sizes = np.bincount(r.ravel())
    profile = counts / bin_sizes
    return profile


def psd2(img):
    return np.abs(np.fft.fftshift(np.fft.fft2(img))) ** 2


radial_psd2 = compose(radial_profile, psd2)


def sharpness(img, q=0.999, radius=1):
    # FROM: https://stackoverflow.com/questions/7765810/is-there-a-way-to-detect-if-an-image-is-blurry/7767755#7767755
    # TODO: vectorize?? normalize?? verify against some kind of ground truth?
    img_blurred = skimage.filters.gaussian(img, radius)
    img_lofg = skimage.filters.laplace(img_blurred)
    return np.percentile(img_lofg, 100 * q)


def hessian_eigenvalues(img, sigma=1.5):
    return hessian_matrix_eigvals(
        hessian_matrix(skimage.img_as_float(img), sigma, mode="nearest", order="rc")
    )


# TODO: use skimage.filters.ridge.compute_hessian_eigenvalues?
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


# FROM: Kovesi, Peter. 2010. Fast Almost-Gaussian Filtering.
# TODO: replace with http://blog.ivank.net/fastest-gaussian-blur.html
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


# FROM: http://emmanuelle.github.io/a-tutorial-on-segmentation.html
def permute_labels(labels):
    label_map = np.concatenate(((0,), np.random.permutation(labels.max()) + 1))
    return label_map[labels]
