import numpy as np
import pandas as pd
import zarr
import skimage
from functools import partial
from trench_detection import get_img_limits
from util import tqdm_auto, map_collections, iterate_get_collection_value


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


def extract_kymograph(img_series, x0, x1):
    num_timepoints = img_series.shape[0]
    xs, ys = coords_along(x0, x1)
    kymo = np.zeros((len(xs), num_timepoints))
    for t in range(num_timepoints):
        kymo[:, t] = img_series[t, ys, xs][::-1]
    return kymo


def bounding_box(points):
    upper_left_x = min(point[0] for point in points)
    upper_left_y = min(point[1] for point in points)
    lower_right_x = max(point[0] for point in points)
    lower_right_y = max(point[1] for point in points)
    return (
        np.array([upper_left_x, upper_left_y]),
        np.array([lower_right_x, lower_right_y]),
    )


def crop_point(x, x_lim, y_lim):
    return np.clip(x, *np.vstack([[x_lim], [y_lim]]).T)


def get_trench_thumbnail(img, trench_points, trench_idx):
    x_lim, y_lim = get_img_limits(img.shape)
    ul, lr = get_trench_bbox(trench_points, trench_idx, x_lim, y_lim)
    return img[ul[1] : lr[1], ul[0] : lr[0]]


def linear_mix(x, y, s):
    return (1 - s) * x + s * y


def get_trench_bbox(trench_points, trench_idx, x_lim, y_lim, overlap=0):
    # separation=0 is the narrowest trench, separation=1 includes neighboring trench midlines
    separation = 0.5 + overlap
    # trench_points[0][trench_idx], trench_points[1][trench_idx]
    num_trenches = min(len(trench_points[0]), len(trench_points[1]))
    if not 0 <= trench_idx < num_trenches:
        raise ValueError("trench index out of bounds")
    points = [trench_points[i][trench_idx] for i in (0, 1)]
    for i in (0, 1):
        if trench_idx == 0:
            x_prev = 2 * trench_points[i][0] - trench_points[i][1]
        else:
            x_prev = trench_points[i][trench_idx - 1]
        x_prev = linear_mix(trench_points[i][trench_idx], x_prev, separation).astype(
            np.int_
        )  # TODO: should we round or truncate? here and below?
        x_prev = crop_point(x_prev, x_lim, y_lim)
        points.append(x_prev)
        if trench_idx + 1 == num_trenches:
            x_next = 2 * trench_points[i][-1] - trench_points[i][-2]
        else:
            x_next = trench_points[i][trench_idx + 1]
        x_next = linear_mix(trench_points[i][trench_idx], x_next, separation).astype(
            np.int_
        )
        x_next = crop_point(x_next, x_lim, y_lim)
        points.append(x_next)
    return bounding_box(points)


def trenchwise_apply(
    img_stack,
    trench_set_points,
    trench_idx,
    func,
    channel_slice=slice(None),
    time_slice=slice(None),
):
    x_lim, y_lim = get_img_limits(img_stack.shape[-2:])
    ul, lr = get_trench_bbox(trench_set_points, trench_idx, x_lim, y_lim)
    if len(img_stack.shape) == 4:
        thumb_stack = (
            img_stack.oindex if isinstance(img_stack, zarr.Array) else img_stack
        )[channel_slice, time_slice, ul[1] : lr[1], ul[0] : lr[0]]
    else:
        thumb_stack = img_stack[..., ul[1] : lr[1], ul[0] : lr[0]]
    res = func(thumb_stack)
    return res


def trenchwise_map(
    img_stack, trench_points, func, progress_bar=tqdm_auto, preload=True, **kwargs
):
    if preload:
        img_stack = (
            img_stack.oindex if isinstance(img_stack, zarr.Array) else img_stack
        )[
            kwargs.get("channel_slice", slice(None)),
            kwargs.get("time_slice", slice(None)),
            :,
            :,
        ]
        del kwargs["channel_slice"]
        del kwargs["time_slice"]
    obj = {
        trench_set_idx: {
            trench_idx: trenchwise_apply(
                img_stack, trench_set_points, trench_idx, func, **kwargs
            )
            for trench_idx in progress_bar(range(len(trench_set_points[0])))
        }
        for trench_set_idx, trench_set_points in progress_bar(trench_points.items())
    }
    representative_obj = iterate_get_collection_value(obj, 2)
    if isinstance(representative_obj, (pd.Series, pd.DataFrame)):
        df = map_collections(partial(pd.concat, axis=1), obj, max_level=2)
        df.columns.set_names("trench_set_idx", level=0, inplace=True)
        df.columns.set_names("trench_idx", level=1, inplace=True)
        return df
    else:
        return obj


def positionwise_trenchwise_map(
    img_group, trench_points_pos, func, positions=None, progress_bar=tqdm_auto, **kwargs
):
    if positions is None:
        positions = trench_points_pos.keys()
    obj = {
        pos: trenchwise_map(
            img_group[pos],
            trench_points_pos[pos],
            func,
            progress_bar=lambda x: x,
            **kwargs,
        )
        for pos in progress_bar(positions)
    }
    if isinstance(iterate_get_collection_value(obj, 1), (pd.Series, pd.DataFrame)):
        df = pd.concat(obj, axis=1)
        df.columns.set_names("position", level=0, inplace=True)
        return df
    else:
        return obj


def image_sharpness(img):
    # FROM: https://stackoverflow.com/questions/7765810/is-there-a-way-to-detect-if-an-image-is-blurry/7767755#7767755
    # TODO: vectorize?? normalize?? verify against some kind of ground truth?
    img_blurred = skimage.filters.gaussian(img, 1)
    img_lofg = skimage.filters.laplace(img_blurred)
    return np.percentile(img_lofg, 99.9)
