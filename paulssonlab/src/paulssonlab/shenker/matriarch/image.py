import numpy as np
import scipy
import skimage.morphology
from util import repeat_apply
import numba
from cytoolz import compose


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


@numba.jit(nopython=True)
def _hough_line_intensity(accumulator, img, theta, diagonal):
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    num_thetas = len(theta)
    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            if img[i, j] == 0:  # optimization for boolean input images
                continue
            for k in range(num_thetas):
                rho = int(np.round_(i * sin_theta[k] + j * cos_theta[k])) + diagonal
                accumulator[rho, k] += img[i, j]


def remove_large_objects(labeled_img, max_size, in_place=False):
    if not in_place:
        labeled_img = labeled_img.copy()
    component_sizes = np.bincount(labeled_img.ravel())
    too_big = component_sizes > max_size
    too_big_mask = too_big[labeled_img]
    labeled_img[too_big_mask] = 0
    return labeled_img


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
        median = np.median(maxes)
    img_labels = repeat_apply(skimage.morphology.dilation, dilation)(img_labels)
    img[img_labels == 0] = 0
    for idx, label in enumerate(label_index):
        mask = img_labels == label
        if weighted:
            img[mask] *= min(maxes[idx] / median, 1)
        else:
            img[mask] /= maxes[idx]
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


radial_psd2 = compose(image.radial_profile, psd2)


def image_sharpness(img):
    # FROM: https://stackoverflow.com/questions/7765810/is-there-a-way-to-detect-if-an-image-is-blurry/7767755#7767755
    # TODO: vectorize?? normalize?? verify against some kind of ground truth?
    img_blurred = skimage.filters.gaussian(img, 1)
    img_lofg = skimage.filters.laplace(img_blurred)
    return np.percentile(img_lofg, 99.9)


def hessian_eigenvalues(img):
    I = skimage.filters.gaussian(img, 1.5)
    I_x = skimage.filters.sobel_h(I)
    I_y = skimage.filters.sobel_v(I)
    I_xx = skimage.filters.sobel_h(I_x)
    I_xy = skimage.filters.sobel_v(I_x)
    I_yx = skimage.filters.sobel_h(I_y)
    I_yy = skimage.filters.sobel_v(I_y)
    kappa_1 = (I_xx + I_yy) / 2
    kappa_2 = (np.sqrt((I_xx + I_yy) ** 2 - 4 * (I_xx * I_yy - I_xy * I_yx))) / 2
    k1 = kappa_1 + kappa_2
    k2 = kappa_1 - kappa_2
    k1[np.isnan(k1)] = 0
    k2[np.isnan(k2)] = 0
    return k1, k2
