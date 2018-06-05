import numpy as np
import numba


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
