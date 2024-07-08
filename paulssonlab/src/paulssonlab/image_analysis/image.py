import warnings
from functools import partial

import numba
import numpy as np
import pandas as pd
import scipy
import scipy.ndimage
import skimage.morphology
import skimage.transform
from cytoolz import compose
from matplotlib.colors import hex2color
from skimage.feature import hessian_matrix, hessian_matrix_eigvals

from paulssonlab.image_analysis.blur import scipy_box_blur
from paulssonlab.image_analysis.geometry import get_image_limits
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


def scale_image(img, scale=False, downsample=None, squeeze=True):
    if squeeze:
        # needed for downsample to work if we have extraneous singleton dimensions up front
        img = img.squeeze()
    if downsample is not None:
        img = img[::downsample, ::downsample, ...]
    if scale is True:
        img_min = np.nanmin(img)
        img = (img - img_min) / (np.nanmax(img) - img_min)
    elif scale:
        img = img - np.nanmin(img)
        img = img / np.nanpercentile(img, scale * 100)
    img = np.clip(img, 0, 1)
    return img


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


def hough_bounds(shape, theta):
    theta = (theta + np.pi / 2) % np.pi - np.pi / 2
    x_lim, y_lim = get_image_limits(shape)
    x_min, x_max = x_lim
    y_min, y_max = y_lim
    if theta < 0:
        rho_min = y_max * np.sin(theta)
        rho_max = x_max * np.cos(theta)
    else:
        rho_min = 0
        rho_max = x_max * np.cos(theta) + y_max * np.sin(theta)
    return rho_min, rho_max


# FROM: https://alyssaq.github.io/2014/understanding-hough-transform/
def hough_line_intensity(img, theta=None):
    if theta is None:
        theta = np.linspace(np.deg2rad(-90), np.deg2rad(90), 180, endpoint=False)
    else:
        theta = np.asarray(theta)
    width, height = img.shape
    diagonal = int(np.ceil(np.sqrt(width**2 + height**2)))
    rho = np.linspace(-diagonal, diagonal, diagonal * 2)
    accumulator = np.zeros((len(rho), len(theta)), dtype=_accumulator_dtype(img.dtype))
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
                rho = int(np.floor(i * sin_theta[k] + j * cos_theta[k])) + diagonal
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


def hessian_eigenvalues(img, sigma=1.5):
    return hessian_matrix_eigvals(
        hessian_matrix(
            skimage.img_as_float(img),
            sigma,
            mode="nearest",
            order="rc",
            use_gaussian_derivatives=True,
        )
    )


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


class TranslationTransform(skimage.transform.EuclideanTransform):
    def estimate(self, src, dst):
        translation = (dst - src).mean(axis=0)
        self.params[0 : self.dimensionality, self.dimensionality] = translation
        return True


# SEE: https://math.stackexchange.com/a/2007723
def _quadratic_root(a, b, c):
    # we've chosen the correct root for for radial_transform (smallest positive root)
    return 2 * c / (-b + np.sqrt(b**2 - 4 * a * c))


# SEE: http://www.cs.ait.ac.th/~mdailey/papers/Bukhari-RadialDistortion.pdf
# Bukhari, F., & Dailey, M. N. (2013). Automatic radial distortion estimation from a single image. Journal of mathematical imaging and vision, 45, 31-45.
def radial_distortion(coords, k1, input_center=None, output_center=None):
    if input_center is not None:
        coords = coords - np.asarray(input_center)
    r_d = np.sqrt((coords**2).sum(axis=-1))[..., np.newaxis]
    new_coords = coords / (1 + k1 * r_d**2)
    if output_center is not None:
        new_coords += np.asarray(output_center)
    return new_coords


# SEE: http://www.cs.ait.ac.th/~mdailey/papers/Bukhari-RadialDistortion.pdf
# Bukhari, F., & Dailey, M. N. (2013). Automatic radial distortion estimation from a single image. Journal of mathematical imaging and vision, 45, 31-45.
def radial_distortion_inverse(coords, k1, input_center=None, output_center=None):
    if output_center is not None:
        coords = coords - np.asarray(output_center)
    r_u = np.sqrt((coords**2).sum(axis=-1))[..., np.newaxis]
    r_d = _quadratic_root(k1 * r_u, -1, r_u)
    with np.errstate(invalid="ignore"):
        # center pixel has r_u = 0, which raises a warning (that we silence)
        new_coords = (r_d / r_u) * coords
    if input_center is not None:
        new_coords += np.asarray(input_center)
    return new_coords


def radial_distortion_output(k1, input_shape, input_center):
    input_shape = np.asarray(input_shape)
    input_center = np.asarray(input_center)
    # we don't use radial_distortion so that coords**2 in denominator
    # is elementwise (no summing to get euclidean distance)
    coords = np.stack((input_center, input_shape[::-1] - input_center))
    output_center, output_lr = coords / (1 + k1 * coords**2)
    output_shape = np.ceil(output_center + output_lr).astype(np.int32)[::-1]
    return output_shape, output_center


def center_from_shape(shape):
    # -1 because skimage interprets whole integer coördinates
    # as the centers of pixels. i.e., (0,0) is the upper-leftmost pixel
    return (np.array(shape)[::-1] - 1) / 2


def radial_distortion_coords(
    k1,
    input_shape,
    input_center=None,
    output_shape=None,
    output_center=None,
):
    if input_center is None:
        input_center = center_from_shape(input_shape)
    else:
        input_center = np.asarray(input_center)
    if (
        num_output_kwargs := sum([output_center is not None, output_shape is not None])
    ) == 1:
        raise ValueError(
            "either both or neither of output_center and output_shape must be specified"
        )
    elif num_output_kwargs == 0:
        output_shape, output_center = radial_distortion_output(
            k1, input_shape, input_center
        )
    return skimage.transform.warp_coords(
        partial(
            radial_distortion_inverse,
            k1=k1,
            input_center=input_center,
            output_center=output_center,
        ),
        output_shape,
    )


def correct_radial_distortion(
    img,
    k1=None,
    coords=None,
    input_center=None,
    output_shape=None,
    output_center=None,
    **kwargs,
):
    if coords is None:
        if k1 is None:
            raise ValueError(
                "either k1 or coords (pre-computed with radial_distortion_coords) must be specified"
            )
        coords = radial_distortion_coords(
            k1,
            img.shape,
            input_center=input_center,
            output_shape=output_shape,
            output_center=output_center,
        )
    return scipy.ndimage.map_coordinates(
        img,
        coords,
        **{**dict(order=1, mode="constant", cval=0, prefilter=False), **kwargs},
    )


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


def unstack(ary, axis=1):
    ary = np.moveaxis(ary, 0, axis - 1)
    return ary.reshape(*ary.shape[: axis - 1], -1, *ary.shape[axis + 1 :])


def pad_and_stack(arys, fill_value=0):
    shape = np.max([ary.shape for ary in arys], axis=0)
    return np.stack(
        [
            np.pad(
                ary,
                ((shape[0] - ary.shape[0], 0), (shape[1] - ary.shape[1], 0)),
                constant_values=fill_value,
            )
            for ary in arys
        ]
    )


def pad_unstack(arys):
    return unstack(pad_and_stack(arys))


def unstack_multichannel(arys, colors=None, scale=True, axis=2):
    imgs = []
    for i in range(arys.shape[0]):
        img = unstack(arys[i], axis=axis).T
        if scale:
            img = scale_image(img, scale=scale)
        imgs.append(img)
    if colors is not None:
        return colorize(imgs, colors, scale=False)
    else:
        return np.concatenate(imgs, axis=-1)


def pad(ary, shape, fill_value=None):
    if fill_value is None:
        if np.issubdtype(ary.dtype, np.integer):
            fill_value = 0
        else:
            fill_value = np.nan
    dtype = np.promote_types(ary.dtype, np.min_scalar_type(fill_value))
    ary = ary.astype(dtype)
    shape = tuple(shape)
    if ary.shape == shape:
        return ary
    else:
        return np.pad(
            ary,
            [(0, max(goal - current, 0)) for goal, current in zip(shape, ary.shape)],
            constant_values=fill_value,
        )


def crop_to_mask(ary, mask=np.isfinite):
    if callable(mask):
        mask = mask(ary)
    for i in range(ary.ndim):
        axis_mask = np.any(mask, axis=tuple([j for j in range(mask.ndim) if j != i]))
        ary = ary[(*[slice(None)] * i, axis_mask)]
    return ary


def colorize(imgs, hexcolors, scale=True):
    colors = [hex2color(hexcolor) for hexcolor in hexcolors]
    return _colorize(imgs, colors, scale=scale)


def _colorize(channel_imgs, colors, scale=True):
    if len(channel_imgs) != len(colors):
        raise ValueError("expecting equal numbers of channels and colors")
    num_channels = len(channel_imgs)
    if scale:
        scaled_imgs = [
            channel_imgs[i] / np.percentile(channel_imgs[i], 99.9)
            for i in range(num_channels)
        ]
        for scaled_img in scaled_imgs:
            np.clip(scaled_img, 0, 1, scaled_img)  # clip in place
    else:
        scaled_imgs = channel_imgs
    imgs_to_combine = [
        scaled_img[:, :, np.newaxis] * np.array(color)
        for scaled_img, color in zip(scaled_imgs, colors)
    ]
    if not len(imgs_to_combine):
        return np.ones(channel_imgs[0].shape)  # white placeholder
    img = imgs_to_combine[0]
    for img2 in imgs_to_combine[1:]:
        img = 1 - (1 - img) * (1 - img2)
    return img


def power_law_composite(imgs, power=0.5, scale=True, dtype=None):
    return np.mean(
        [scale_image(img, scale=scale, squeeze=False) ** power for img in imgs],
        axis=0,
        dtype=dtype,
    )
