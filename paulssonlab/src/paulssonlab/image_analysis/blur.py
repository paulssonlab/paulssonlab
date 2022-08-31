import numpy as np
import scipy.ndimage as ndi
import numba

# SEE:
# - https://stackoverflow.com/questions/98359/fastest-gaussian-blur-implementation
# - https://github.com/bfraboni/FastGaussianBlur
# - http://blog.ivank.net/fastest-gaussian-blur.html
# - https://dsp.stackexchange.com/questions/50576/fastest-available-algorithm-to-blur-an-image-low-pass-filter
# - https://github.com/pelson/antigrain/blob/64c9125e2b350a422c08d7fa8fff023400ad3f9f/agg-2.4/include/agg_blur.h


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
def scipy_box_blur_2d(ary, sigma, n=3, mode="nearest"):
    m, w_l, w_u = box_blur_params(sigma, n)
    for i in range(m):
        ary = ndi.filters.uniform_filter(ary, w_l, mode=mode)
    for i in range(n - m):
        ary = ndi.filters.uniform_filter(ary, w_u, mode=mode)
    return ary


def _box_blur(filter_1d, ary, sigma, n, mode="nearest"):
    m, w_l, w_u = box_blur_params(sigma, n)
    for i in range(m):
        ary = filter_1d(ary, w_l, mode=mode)
    for i in range(n - m):
        ary = filter_1d(ary, w_u, mode=mode)
    # SEE https://stackoverflow.com/a/62523103 for how to do transpose in place
    # this would save memory, but is 10x slower
    # see also fastremap.transpose (https://github.com/seung-lab/fastremap)
    # turns out that it's faster to avoid copy (physical transpose) here and below
    ary = ary.transpose()
    for i in range(m):
        ary = filter_1d(ary, w_l, mode=mode)
    for i in range(n - m):
        ary = filter_1d(ary, w_u, mode=mode)
    ary = ary.transpose()
    return ary


# same as scipy_box_blur_2d but does all horizontal passes first, then all vertical passes
def scipy_box_blur(ary, sigma, n=3, mode="nearest"):
    return _box_blur(ndi.filters.uniform_filter1d, ary, sigma, n, mode=mode)


def _numba_uniform_filter1d(ary, w, mode="nearest"):
    # TODO: implement mode
    # copy one row at a time into original array (in place operation)
    normalization = 1 / (2 * r + 1)
    for i in range(ary.shape[-2]):
        for j in range(w):
            total += ary[i, j]
        for j in range(w + 1):
            val += ary[i, j] - fv


# FROM: http://blog.ivank.net/fastest-gaussian-blur.html
def numba_box_blur(ary, sigma, n=3, mode="nearest"):
    return _box_blur(_numba_uniform_filter1d, sigma, n, mode=mode)


@numba.njit(fastmath=True, error_model="numpy")
def _numba_stack_blur1d(ary, radius):
    dtype = ary.dtype
    width = ary.shape[-1]
    width_minus1 = width - 1
    height = ary.shape[-2]
    radius_plus1 = radius + 1
    div = 2 * radius + 1
    div_sum = (radius + 1) * (radius + 1)
    sum_factor = (radius + 1) * (radius + 2) // 2
    new_row = np.empty(width, dtype=dtype)
    stack = np.empty(div, dtype=dtype)
    for y in range(height):
        sum_ = 0
        sum_in = 0
        sum_out = 0
        value = ary[y, 0]
        sum_out = radius_plus1 * value
        sum_ += sum_factor * value
        for i in range(radius_plus1):
            stack[i] = value
        for i in range(1, radius_plus1):
            value = ary[y, width_minus1 if i > width_minus1 else i]
            stack[i + radius] = value
            sum_ += value * (radius_plus1 - i)
            sum_in += value
        stack_idx = radius
        for x in range(width):
            new_row[x] = np.round(sum_ / div_sum)
            sum_ -= sum_out
            stack_start = (stack_idx + div - radius) % div
            stack_value = stack[stack_start]
            sum_out -= stack_value
            xp = (x + radius_plus1) % width_minus1
            value = ary[y, xp]
            stack[stack_start] = value
            sum_in += value
            sum_ += sum_in
            stack_idx = (stack_idx + 1) % div
            stack_value = stack[stack_idx]
            sum_out += stack_value
            sum_in -= stack_value
        # copy new_row to row
        for x in range(width):
            ary[y, x] = new_row[x]


def numba_stack_blur(ary, sigma, mode="nearest"):
    radius = sigma  # TODO: what radius gives best approximation of gaussian blur?
    if mode != "nearest":
        raise ValueError("only mode 'nearest' is supported")
    _numba_stack_blur1d(ary, radius)
    ary = ary.transpose()
    _numba_stack_blur1d(ary, radius)
    ary = ary.transpose()
    return ary
