import numpy as np


def get_image_limits(shape):
    x_min = y_min = 0
    y_max, x_max = shape  # note order
    # TODO: what convention should we use, should max be inclusive??
    # x_max = img.shape[0] - 1
    # y_max = img.shape[1] - 1
    x_lim = (x_min, x_max)
    y_lim = (y_min, y_max)
    return x_lim, y_lim


# TODO: obsolete
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
    x_lim, y_lim = get_image_limits(img.shape)
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
