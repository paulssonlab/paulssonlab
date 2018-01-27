import numpy as np
import scipy
import skimage
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
from utils import get_if_not_none


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
    return label_img


def drop_rare_labels(labels):
    counter = Counter(labels)
    total = sum(counter)
    good_labels = []
    for label, count in counter.iteritems():
        print(count / total)
        if count / total > 0.01:
            good_labels.append(label)
    return good_labels


def get_img_limits(img):
    x_min = y_min = 0
    x_max, y_max = img.shape
    # TODO: what convention should we use, should max be inclusive??
    # x_max = img.shape[0] - 1
    # y_max = img.shape[1] - 1
    x_lim = (x_min, x_max)
    y_lim = (y_min, y_max)
    return x_lim, y_lim


def detect_rotation(bin_img, diagnostics=None):
    h, theta, d = skimage.transform.hough_line(bin_img)
    abs_diff_h = np.diff(h.astype(np.int32), axis=1).var(axis=0)
    theta_idx = abs_diff_h.argmax()
    angle1 = theta[theta_idx]
    h2, theta2, d2 = skimage.transform.hough_line(
        bin_img, theta=np.linspace(0.9 * angle1, 1.1 * angle1, 200)
    )
    abs_diff_h2 = np.diff(h2.astype(np.int32), axis=1).var(axis=0)
    theta_idx2 = abs_diff_h2.argmax()
    angle2 = theta2[theta_idx2]
    d_profile = h2[:, theta_idx2].astype(np.int32)
    freqs = np.abs(np.fft.fft(d_profile))
    peak_idxs = peakutils.indexes(d_profile, thres=0.4, min_dist=5)
    peaks = d2[peak_idxs]
    spacing = scipy.stats.mode(np.diff(peaks)).mode[0]
    return np.pi / 2 - angle2, peaks


def get_rough_spacing(dists):
    spacing = scipy.stats.mode(np.diff(dists).astype(int)).mode[0]
    return spacing


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


def get_anchors(theta, x_lim, y_lim):
    x_min = np.array([x_lim[0], y_lim[0]])
    x_max = np.array([x_lim[1], y_lim[1]])
    x0 = x_min + (x_max - x_min) / 2
    anchor0 = edge_point(x0, theta, x_lim, y_lim)
    anchor1 = edge_point(x0, theta + np.pi, x_lim, y_lim)
    return anchor0, anchor1


def detect_trench_region(bin_img, theta):
    x_lim, y_lim = get_img_limits(bin_img)
    anchor0, anchor1 = get_anchors(theta, x_lim, y_lim)
    cross_sections = []
    anchors = list(point_linspace(anchor0, anchor1, 40))[3:-3]  # TODO: parameterize
    lines = list(
        line_array(anchors, np.pi / 2 + theta, x_lim, y_lim, bidirectional=True)
    )
    for x0, x1, x2 in lines:
        xs, ys = coords_along(x1, x2)
        cross_sections.append(bin_img[ys, xs])
    cross_section_vars = np.array([cs.var() for cs in cross_sections])
    idx = cross_section_vars.argmax()
    return anchors[idx]


# FROM: https://stackoverflow.com/questions/23815327/numpy-one-liner-for-combining-unequal-length-np-array-to-a-matrixor-2d-array
def stack_jagged(arys, fill=0):
    return np.array(list(zip_longest(*arys, fillvalue=fill))).T


def time_series_periodogram(xs, period_min, period_max, bins=100, diagnostics=None):
    periods = np.linspace(period_min, period_max, bins)
    # std = ((dxs + periods[:,np.newaxis]/2) % periods[:,np.newaxis]).std(axis=1)
    std = scipy.stats.iqr(((xs) % periods[:, np.newaxis]), axis=1) / periods
    period_idx = std.argmin()
    period = periods[period_idx]
    if diagnostics is not None:
        highlighted_point = hv.Points([(period, std[period_idx])]).opts(
            style={"Points": {"size": 10, "color": "red"}}
        )
        diagnostics["periodogram"] = hv.Curve((periods, std)) * highlighted_point
    # plt.figure()
    # plt.plot(periods, std)
    # plt.scatter([period], [std[period_idx]], c='r')
    return period


def detect_periodic_peaks(signal, min_dist=10, max_dist=50, diagnostics=None):
    idxs = peakutils.indexes(signal, thres=0.2, min_dist=min_dist)
    # xs = peakutils.interpolate(np.arange(len(signal)), signal, ind=idxs)
    xs = idxs
    # plt.figure(figsize=(16,8))
    # plt.plot()
    # plt.plot(np.arange(len(signal)), signal)
    # plt.scatter(xs, signal[xs.astype(np.int_)], c='r')
    dxs = np.diff(xs)
    # period_min = np.percentile(dxs, 10)
    # period_max = dxs.max()
    period_min = min_dist
    period_max = max_dist
    num_periods = 100
    period = time_series_periodogram(
        xs,
        period_min,
        period_max,
        bins=num_periods,
        diagnostics=get_if_not_none(diagnostics, "periodogram_1"),
    )
    # period2 = period
    period2 = time_series_periodogram(
        xs,
        period * 0.98,
        period * 1.02,
        num_periods,
        diagnostics=get_if_not_none(diagnostics, "periodogram_2"),
    )
    offsets = np.linspace(0, period2, num_periods)
    offset_idxs = (
        np.arange(0, len(signal) - period2, period2) + offsets[:, np.newaxis]
    ).astype(np.int_)
    objective = signal[offset_idxs].sum(axis=1)
    offset_idx = objective.argmax()
    offset = offsets[offset_idx]
    # plt.figure()
    # plt.plot(offsets, objective)
    # plt.scatter([offset], [objective[offset_idx]], c='r')
    return period2, offset


def detect_trench_anchors(img, t0, theta, diagnostics=None):
    x_lim, y_lim = get_img_limits(img)
    x1 = edge_point(t0, theta - np.pi / 2, x_lim, y_lim)
    x2 = edge_point(t0, theta + np.pi / 2, x_lim, y_lim)
    xs, ys = coords_along(x1, x2)
    profile = img[ys, xs]
    period, offset = detect_periodic_peaks(profile, diagnostics=diagnostics)
    idxs = np.arange(offset, len(profile), period).astype(np.int_)
    # plt.figure(figsize=(16,8))
    # plt.plot(profile)
    # plt.scatter(idxs, profile[idxs], c='r')
    return np.vstack((xs[idxs], ys[idxs])).T


def _detect_trench_end(img, anchors, theta, diagnostics=None):
    x_lim, y_lim = get_img_limits(img)
    xss = []
    yss = []
    trench_profiles = []
    for anchor in anchors:
        x_end = edge_point(anchor, theta, x_lim, y_lim)
        xs, ys = coords_along(anchor, x_end)
        xss.append(xs)
        yss.append(ys)
        trench_profiles.append(img[ys, xs])
    # plt.figure(figsize=(8,8))
    # for trench_profile in trench_profiles:
    #     plt.plot(trench_profile)
    stacked_profile = np.percentile(stack_jagged(trench_profiles), 80, axis=0)
    # cum_profile = np.cumsum(stacked_profile)
    # cum_profile /= cum_profile[-1]
    # end = np.where(cum_profile > 0.8)[0][0]
    stacked_profile_diff = holo_diff(1, stacked_profile)
    end = stacked_profile_diff.argmin()
    # plt.figure(figsize=(8,8))
    # plt.plot(stacked_profile)
    # plt.axvline(end, c='r')
    # plt.figure(figsize=(8,8))
    # plt.plot(stacked_profile_diff, color='g')
    # plt.axvline(end, c='r')
    end_points = []
    for xs, ys in zip(xss, yss):
        idx = end
        if len(xs) <= end:
            idx = -1
        end_points.append((xs[idx], ys[idx]))
    return np.array(end_points)


def detect_trench_ends(img, bin_img, anchors, theta, diagnostics=None):
    img_masked = np.where(
        skimage.morphology.binary_dilation(bin_img), img, np.percentile(img, 5)
    )
    top_points = _detect_trench_end(
        img_masked, anchors, theta, diagnostics=get_if_not_none(diagnostics, "top")
    )
    bottom_points = _detect_trench_end(
        img_masked,
        anchors,
        theta + np.pi,
        diagnostics=get_if_not_none(diagnostics, "bottom"),
    )
    # plt.figure(figsize=(12,12))
    # plt.imshow(img_masked)
    # plt.scatter(*anchors.T, s=3, c='w')
    # plt.scatter(*top_points.T, s=3, c='g')
    # plt.scatter(*bottom_points.T, s=3, c='r')
    return top_points, bottom_points


def detect_trenches(img, bin_img, theta, diagnostics=None):
    t0 = detect_trench_region(bin_img, theta)
    trench_anchors = detect_trench_anchors(img, t0, theta, diagnostics=diagnostics)
    trench_points = detect_trench_ends(
        img, bin_img, trench_anchors, theta, diagnostics=diagnostics
    )
    return trench_points


def _label_for_trenches(img_series, channel):
    # img = img_series[::10].max(axis=0)
    img = img_series[channel, 30]  # TODO
    # TODO: need rotation-invariant detrending
    img = img - np.percentile(img, 3, axis=1)[:, np.newaxis]
    img_thresh = img > skimage.filters.threshold_otsu(img)
    img_labels = label_binary_image(img_thresh)
    return img, img_labels


def get_trenches(img_series, channel, diagnostics=None):
    img, img_labels = _label_for_trenches(img_series, channel)
    if diagnostics is not None:
        diagnostics["labeled_image"] = datashader.regrid(hv.Image(img))
    max_label = img_labels.max()
    trenches = {}
    for label in range(1, max_label + 1):  # TODO: this relies on background == 0
        trenches[label] = _get_trench_set(
            img,
            img_labels == label,
            diagnostics=get_if_not_none(diagnostics, "label_{}".format(label)),
        )
    return trenches


def _get_trench_set(img, img_mask, diagnostics=None):
    theta, dists = detect_rotation(img_mask, diagnostics=diagnostics)
    trench_points = detect_trenches(img, img_mask, theta, diagnostics=diagnostics)
    return trench_points
