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


def detect_trench_region(
    bin_img, theta, margin=3, num_region_slices=100, smooth=2, diagnostics=None
):
    x_lim, y_lim = get_image_limits(bin_img.shape)
    anchor0, anchor1 = get_edge_points(theta, x_lim, y_lim)
    cross_sections = []
    anchors = list(point_linspace(anchor0, anchor1, num_region_slices))[margin:-margin]
    lines = list(
        line_array(anchors, np.pi / 2 + theta, x_lim, y_lim, bidirectional=True)
    )
    for x0, x1, x2 in lines:
        xs, ys = coords_along(x1, x2)
        cross_sections.append(bin_img[ys, xs])
    if diagnostics is not None:
        # anchors = np.vstack((xs[idxs], ys[idxs])).T
        # for l in lines:
        #    print(l, (l[1][0], l[2][0]), (l[1][1], l[2][1]))
        region_slices_plot = hv.Overlay.from_values(
            [hv.Path(((l[1][0], l[2][0]), (l[1][1], l[2][1]))) for l in lines]
        )
        diagnostics["region_slices"] = RevImage(bin_img) * region_slices_plot
        diagnostics["cross_sections"] = hv.Overlay.from_values(
            [hv.Curve(cs) for cs in cross_sections]
        )
    cross_section_vars = np.array([cs.var() for cs in cross_sections])
    if smooth:
        cross_section_vars_smoothed = scipy.ndimage.filters.gaussian_filter1d(
            cross_section_vars, smooth
        )
    else:
        cross_section_vars_smoothed = cross_section_vars
    idx = cross_section_vars_smoothed.argmax()
    if diagnostics is not None:
        cross_section_vars_plot = hv.Curve(cross_section_vars)
        if smooth:
            cross_section_vars_plot *= hv.Curve(cross_section_vars_smoothed).options(
                color="cyan"
            )
        cross_section_vars_plot *= hv.VLine(idx).options(color="red")
        diagnostics["cross_section_vars"] = cross_section_vars_plot
    return anchors[idx]


# FROM: https://stackoverflow.com/questions/23815327/numpy-one-liner-for-combining-unequal-length-np-array-to-a-matrixor-2d-array
def stack_jagged(arys, fill=0):
    return np.array(list(zip_longest(*arys, fillvalue=fill))).T


# TODO: unused
def discrete_periodogram(xs, period_min, period_max, bins=1000, diagnostics=None):
    periods = np.linspace(period_min, period_max, bins)
    std = scipy.stats.iqr(((xs) % periods[:, np.newaxis]), axis=1) / periods
    period_idx = std.argmin()
    period = periods[period_idx]
    if diagnostics is not None:
        diagnostics["search_range"] = (period_min, period_max)
        diagnostics["bins"] = bins
        diagnostics["period"] = period
        highlighted_point = hv.Scatter([(period, std[period_idx])]).options(
            size=5, color="red"
        )
        diagnostics["periodogram"] = hv.Curve((periods, std)) * highlighted_point
    return period


# TODO: unused
def detect_periodic_peaks(
    signal, threshold=0.2, min_dist=5, period_min=None, max_dist=None, diagnostics=None
):
    if diagnostics is not None:
        # plot signal without highlighted points in case peakutils hangs
        diagnostics["peaks"] = hv.Curve((range(len(signal)), signal))
    if np.all(signal == signal[0]):
        raise ValueError("attempting to detect periodic peaks in a constant signal")
    idxs = peakutils.indexes(signal, thres=threshold, min_dist=min_dist)
    # xs = peakutils.interpolate(np.arange(len(signal)), signal, ind=idxs)
    xs = idxs
    if diagnostics is not None:
        highlighted_points = hv.Scatter((xs, signal[xs.astype(np.int_)])).options(
            size=5, color="red"
        )
        diagnostics["peaks"] = (
            hv.Curve((range(len(signal)), signal)) * highlighted_points
        )
    dxs = np.diff(xs)
    if min_dist is None:
        min_dist = max(dxs.min() - 1, 2)
    if max_dist is None:
        max_dist = dxs.max() + 1
    num_periods = 1000
    period = discrete_periodogram(
        xs,
        min_dist,
        max_dist,
        bins=num_periods,
        diagnostics=getitem_if_not_none(diagnostics, "periodogram_1"),
    )
    period2 = discrete_periodogram(
        xs,
        period * 0.98,
        period * 1.02,
        num_periods,
        diagnostics=getitem_if_not_none(diagnostics, "periodogram_2"),
    )
    offsets = np.linspace(0, period2, num_periods)
    offset_idxs = (
        np.arange(0, len(signal) - period2, period2) + offsets[:, np.newaxis]
    ).astype(np.int_)
    objective = signal[offset_idxs].sum(axis=1)
    offset_idx = objective.argmax()
    offset = offsets[offset_idx]
    if diagnostics is not None:
        highlighted_point = hv.Scatter([(offset, objective[offset_idx])]).options(
            size=5, color="red"
        )
        diagnostics["offsets"] = hv.Curve((offsets, objective)) * highlighted_point
    return period2, offset


def detect_trench_anchors(img, t0, theta, diagnostics=None):
    x_lim, y_lim = get_image_limits(img.shape)
    x1 = edge_point(t0, theta - np.pi / 2, x_lim, y_lim)
    x2 = edge_point(t0, theta + np.pi / 2, x_lim, y_lim)
    xs, ys = coords_along(x1, x2)
    profile = img[ys, xs]  # .astype(np.int_)
    if diagnostics is not None:
        diagnostics["trench_anchor_profile"] = hv.Curve(profile)
        diagnostics["trench_slice"] = RevImage(img) * hv.Path(
            ((x1[0], x2[0]), (x1[1], x2[1]))
        )
    period, offset = detect_periodic_peaks(profile, diagnostics=diagnostics)
    idxs = np.arange(offset, len(profile), period).astype(np.int_)
    anchors = np.vstack((xs[idxs], ys[idxs])).T
    if diagnostics is not None:
        highlighted_points = hv.Scatter((idxs, profile[idxs])).options(
            size=5, color="red"
        )
        diagnostics["trench_anchor_profile"] = hv.Curve(profile) * highlighted_points
        anchor_points_plot = hv.Points(anchors).options(size=3, color="white")
        diagnostics["anchor_points"] = RevImage(img) * anchor_points_plot
    return anchors


def _detect_trench_end(img, anchors, theta, margin=15, threshold=0.8, diagnostics=None):
    x_lim, y_lim = get_image_limits(img.shape)
    xss = []
    yss = []
    trench_profiles = []
    for anchor in anchors:
        x_end = edge_point(anchor, theta, x_lim, y_lim)
        xs, ys = coords_along(anchor, x_end)
        xss.append(xs)
        yss.append(ys)
        trench_profiles.append(img[ys, xs])
    if diagnostics is not None:
        diagnostics["trench_profiles"] = hv.Overlay.from_values(
            [hv.Curve(tp) for tp in trench_profiles]
        )
        diagnostics["threshold"] = threshold
    stacked_profile = np.percentile(
        stack_jagged(trench_profiles), threshold * 100, axis=0
    )
    # cum_profile = np.cumsum(stacked_profile)
    # cum_profile /= cum_profile[-1]
    # end = np.where(cum_profile > 0.8)[0][0]
    stacked_profile_diff = holo_diff(1, stacked_profile)
    end = stacked_profile_diff.argmin() + margin
    if diagnostics is not None:
        diagnostics["margin"] = margin
        diagnostics["stacked_profile"] = hv.Curve(stacked_profile) * hv.VLine(
            end
        ).options(color="red")
        diagnostics["stacked_profile_diff"] = hv.Curve(stacked_profile_diff) * hv.VLine(
            end
        ).options(color="green")
    end_points = []
    for xs, ys in zip(xss, yss):
        idx = end
        if len(xs) <= end:
            idx = -1
        end_points.append((xs[idx], ys[idx]))
    return np.array(end_points)


def detect_trench_ends(img, anchors, theta, diagnostics=None):
    top_points = _detect_trench_end(
        img, anchors, theta, diagnostics=getitem_if_not_none(diagnostics, "top")
    )
    bottom_points = _detect_trench_end(
        img,
        anchors,
        theta + np.pi,
        diagnostics=getitem_if_not_none(diagnostics, "bottom"),
    )
    if diagnostics is not None:
        anchor_points_plot = hv.Points(anchors).options(size=3, color="white")
        top_points_plot = hv.Points(top_points).options(size=3, color="green")
        bottom_points_plot = hv.Points(bottom_points).options(size=3, color="red")
        # diagnostics['image_with_trenches'] = datashader.regrid(RevImage(img_masked), aggregator='first').redim.range(z=(0,img_masked.max())) * anchor_points_plot * top_points_plot * bottom_points_plot
        diagnostics["image_with_trenches"] = (
            RevImage(img) * anchor_points_plot * top_points_plot * bottom_points_plot
        )
    return top_points, bottom_points
