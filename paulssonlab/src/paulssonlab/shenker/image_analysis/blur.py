import numpy as np
import scipy.ndimage as ndi
import numba

# SEE:
# - https://stackoverflow.com/questions/98359/fastest-gaussian-blur-implementation
# - https://github.com/bfraboni/FastGaussianBlur
# - http://blog.ivank.net/fastest-gaussian-blur.html
# - https://dsp.stackexchange.com/questions/50576/fastest-available-algorithm-to-blur-an-image-low-pass-filter


def box_blur_params(sigma, n):
    w_ideal = np.sqrt((12 * sigma**2 / n) + 1)
    w_l = int(np.floor(w_ideal))
    if w_l % 2 == 0:
        w_l -= 1
    w_u = w_l + 2
    m = (12 * sigma**2 - n * w_l**2 - 4 * n * w_l - 3 * n) / (-4 * w_l - 4)
    m = round(m)
    return m, w_l, w_u


# FROM: Kovesi, Peter. 2010. Fast Almost-Gaussian Filtering.
def scipy_box(ary, sigma, n=3, mode="nearest"):
    m, w_l, w_u = box_blur_params(sigma, n)
    for i in range(m):
        ary = ndi.filters.uniform_filter(ary, w_l, mode=mode)
    for i in range(n - m):
        ary = ndi.filters.uniform_filter(ary, w_u, mode=mode)
    return ary


def scipy_box_transpose(ary, sigma, n=3, mode="nearest"):
    m, w_l, w_u = box_blur_params(sigma, n)
    for i in range(m):
        ary = ndi.filters.uniform_filter1d(ary, w_l, mode=mode)
    for i in range(n - m):
        ary = ndi.filters.uniform_filter1d(ary, w_u, mode=mode)
    # SEE https://stackoverflow.com/a/62523103 for how to do transpose in place
    # this would save memory, but is 10x slower
    # see also fastremap.transpose (https://github.com/seung-lab/fastremap)
    ary = ary.transpose().copy()
    for i in range(m):
        ary = ndi.filters.uniform_filter1d(ary, w_l, mode=mode)
    for i in range(n - m):
        ary = ndi.filters.uniform_filter1d(ary, w_u, mode=mode)
    # TODO: depending on downstream processing, we may or not want to copy here
    ary = ary.transpose().copy()
    return ary


# TODO: numba box blur, stackblur implementations
