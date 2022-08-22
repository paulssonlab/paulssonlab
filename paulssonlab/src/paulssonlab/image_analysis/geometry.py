import numpy as np


def get_image_limits(shape):
    x_min = y_min = 0
    y_max, x_max = shape  # note order
    # TODO: what convention should we use, should max be inclusive??
    x_max -= 1
    y_max -= 1
    x_lim = (x_min, x_max)
    y_lim = (y_min, y_max)
    return x_lim, y_lim


def bounding_box(points):
    upper_left_x = min(point[0] for point in points)
    upper_left_y = min(point[1] for point in points)
    lower_right_x = max(point[0] for point in points)
    lower_right_y = max(point[1] for point in points)
    # return (np.array([upper_left_x, upper_left_y]), np.array([lower_right_x, lower_right_y]))
    return np.array([[upper_left_x, upper_left_y], [lower_right_x, lower_right_y]])


def crop_point(x, x_lim, y_lim):
    return np.clip(x, *np.vstack([[x_lim], [y_lim]]).T)


def get_trench_thumbnail(img, trench_points, trench_idx):
    x_lim, y_lim = get_image_limits(img.shape)
    ul, lr = get_trench_bbox(trench_points, trench_idx, x_lim, y_lim)
    return img[ul[1] : lr[1], ul[0] : lr[0]]


def linear_mix(x, y, s):
    return (1 - s) * x + s * y


def get_trench_bbox(top, bottom, width, trench_idx, x_lim, y_lim):
    half_width = int(np.ceil(width / 2))
    offset = np.array([half_width, 0])
    return np.vstack((top[trench_idx] - offset, bottom[trench_idx] + offset))


# TODO: first and last trenches are broken
def get_periodic_trench_bbox(top, bottom, trench_idx, x_lim, y_lim, overlap=0):
    trench_points = (top, bottom)
    # separation=0 is the narrowest trench, separation=1 includes neighboring trench midlines
    separation = 0.5 + overlap
    if len(top) != len(bottom):
        raise ValueError("top and bottom should have equal lengths")
    num_trenches = len(top)
    if num_trenches <= 1:
        raise ValueError("need at least two trenches to get bboxes")
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
        # if trench_idx == 60:
        #     print("i", i, x_prev)
        points.append(x_prev)
        if trench_idx + 1 == num_trenches:
            x_next = 2 * trench_points[i][-1] - trench_points[i][-2]
        else:
            x_next = trench_points[i][trench_idx + 1]
        x_next = linear_mix(trench_points[i][trench_idx], x_next, separation).astype(
            np.int_
        )
        x_next = crop_point(x_next, x_lim, y_lim)
        # if trench_idx == 60:
        #     print("i", i, x_next)
        points.append(x_next)
    # if trench_idx == 60:
    #     print("bbox", bounding_box(points))
    #     0 / 0
    return bounding_box(points)
