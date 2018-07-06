import numpy as np
import scipy
import skimage
import skimage.morphology
import skimage.segmentation
import sklearn
import sklearn.cluster
from sklearn.preprocessing import StandardScaler
import peakutils
from holoborodko_diff import holo_diff
from collections import defaultdict
from itertools import zip_longest
import matplotlib.pyplot as plt
import holoviews as hv
import holoviews.operation.datashader as datashader
from util import getitem_if_not_none
from ui import RevImage
from image import (
    hough_line_intensity,
    remove_large_objects,
    normalize_componentwise,
    gaussian_box_approximation,
)
from geometry import get_image_limits
import common


def _standardize_cluster_labels(X, fit):
    mean = defaultdict(lambda: 0)
    count = defaultdict(lambda: 0)
    for i in range(len(fit.labels_)):
        mean[fit.labels_[i]] += X[i]
        count[fit.labels_[i]] += 1
    for k, v in mean.items():
        mean[k] = v / count[k]
    label_mapping = dict(zip(mean.keys(), np.lexsort(list(mean.values()))))
    for i, old_label in enumerate(fit.labels_):
        fit.labels_[i] = label_mapping[old_label]


def cluster_binary_image(bin_img):
    X = np.array(np.where(bin_img)).T
    X2 = StandardScaler().fit_transform(X.astype(np.float32))
    fit = sklearn.cluster.MiniBatchKMeans(
        init="k-means++", n_clusters=2, n_init=10, max_no_improvement=10, verbose=0
    )
    fit.fit(X2)
    _standardize_cluster_labels(X, fit)
    return X, fit


def label_binary_image(bin_img):
    X, fit = cluster_binary_image(bin_img)
    label_img = np.zeros_like(bin_img, dtype=np.int8)  # TODO: fixing dtype
    for i in range(len(fit.labels_)):
        label_img[X[i, 0], X[i, 1]] = fit.labels_[i] + 1
    return label_img, np.sort(np.unique(fit.labels_)) + 1


def drop_rare_labels(labels):
    counter = Counter(labels)
    total = sum(counter)
    good_labels = []
    for label, count in counter.iteritems():
        print(count / total)
        if count / total > 0.01:
            good_labels.append(label)
    return good_labels


def point_linspace(anchor0, anchor1, num_points):
    for s in np.linspace(0, 1, num_points)[1:-1]:
        anchor = (1 - s) * anchor0 + s * anchor1
        yield anchor


def coords_along(x0, x1):
    length = int(np.sqrt(np.sum((x1 - x0) ** 2)))
    xs = np.linspace(x0[0], x1[0], length).astype(np.int_)[1:-1]
    ys = np.linspace(x0[1], x1[1], length).astype(np.int_)[1:-1]
    return xs, ys


def edge_point(x0, theta, x_lim, y_lim):
    x_min, x_max = x_lim
    y_min, y_max = y_lim
    theta = theta % (2 * np.pi)
    if 0 <= theta < np.pi / 2:
        corner_x, corner_y = x_min, y_max
    elif np.pi / 2 <= theta < np.pi:
        corner_x, corner_y = x_max, y_max
    elif np.pi <= theta < 3 / 2 * np.pi:
        corner_x, corner_y = x_max, y_min
    elif 3 / 2 * np.pi <= theta <= 2 * np.pi:
        corner_x, corner_y = x_min, y_min
    angle_to_corner = np.arctan2(corner_y - x0[1], x0[0] - corner_x) % (2 * np.pi)
    if (
        (theta >= angle_to_corner and 0 <= theta < np.pi / 2)
        or (theta < angle_to_corner and np.pi / 2 <= theta < np.pi)
        or (theta >= angle_to_corner and np.pi <= theta < 3 / 2 * np.pi)
        or (theta < angle_to_corner and 3 / 2 * np.pi <= theta < 2 * np.pi)
    ):
        # top/bottom
        x1 = np.array([x0[0] - (corner_y - x0[1]) / np.tan(theta), corner_y])
    else:
        # left/right
        x1 = np.array([corner_x, x0[1] - (corner_x - x0[0]) * np.tan(theta)])
    return x1


def line_array(
    anchors, theta, x_lim, y_lim, start=None, stop=None, bidirectional=False
):
    if bidirectional:
        line_array1 = line_array(
            anchors, theta, x_lim, y_lim, start=start, stop=stop, bidirectional=False
        )
        line_array2 = line_array(
            anchors,
            theta + np.pi,
            x_lim,
            y_lim,
            start=start,
            stop=stop,
            bidirectional=False,
        )
        for (x0, x1), (y0, y1) in zip(line_array1, line_array2):
            yield x0, x1, y1
        return
    if start is None:
        start = 0
    if stop is None:
        stop = 0
    if not stop >= start >= 0:
        raise ValueError("need stop >= start >= 0")
    theta = theta % (2 * np.pi)
    for anchor in anchors:
        x0 = anchor
        x1 = edge_point(x0, theta, x_lim, y_lim)
        max_length = np.sqrt(((x1 - x0) ** 2).sum())
        y0, y1 = x0, x1
        if start:
            y0 = min(start / max_length, 1) * (x1 - x0) + x0
        if stop:
            y1 = min(stop / max_length, 1) * (x1 - x0) + x0
        if not np.array_equal(y0, y1):
            yield y0, y1


def get_edge_points(theta, x_lim, y_lim):
    x_min = np.array([x_lim[0], y_lim[0]])
    x_max = np.array([x_lim[1], y_lim[1]])
    x0 = x_min + (x_max - x_min) / 2
    anchor0 = edge_point(x0, theta, x_lim, y_lim)
    anchor1 = edge_point(x0, theta + np.pi, x_lim, y_lim)
    return anchor0, anchor1


# FROM: https://stackoverflow.com/questions/23815327/numpy-one-liner-for-combining-unequal-length-np-array-to-a-matrixor-2d-array
def stack_jagged(arys, fill=0):
    return np.array(list(zip_longest(*arys, fillvalue=fill))).T
