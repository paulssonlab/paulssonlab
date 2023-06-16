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
    return np.array([[upper_left_x, upper_left_y], [lower_right_x, lower_right_y]])


def iter_crops(img, trenches, corner=False):
    index = trenches.index.values
    ul_x = trenches["ul_x"].values
    ul_y = trenches["ul_y"].values
    lr_x = trenches["lr_x"].values
    lr_y = trenches["lr_y"].values
    for i in range(len(index)):
        crop = img[ul_y[i] : lr_y[i] + 1, ul_x[i] : lr_x[i] + 1]
        if corner:
            yield i, crop, np.array([ul_x[i], ul_y[i]])
        else:
            yield i, crop
