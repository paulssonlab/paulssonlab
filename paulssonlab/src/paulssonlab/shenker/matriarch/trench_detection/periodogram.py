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
from itertools import zip_longest
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
from .core import label_binary_image


def find_hough_angle(
    img, theta=None, smooth=4, hough_func=hough_line_intensity, diagnostics=None
):
    if theta is None:
        theta = np.linspace(np.deg2rad(-45), np.deg2rad(45), 90)
    h, theta, d = hough_func(img, theta=theta)
    diff_h = np.diff(h.astype(np.int_), axis=1)
    diff_h_std = diff_h.std(axis=0)  # / diff_h.max(axis=0)
    if smooth:
        diff_h_std_smoothed = scipy.ndimage.filters.gaussian_filter1d(
            diff_h_std, smooth
        )
    else:
        diff_h_std_smoothed = diff_h_std
    theta_idx = diff_h_std_smoothed.argmax()
    angle = theta[theta_idx]
    if diagnostics is not None:
        diagnostics["angle_range"] = (np.rad2deg(theta[0]), np.rad2deg(theta[-1]))
        bounds = (np.rad2deg(theta[0]), d[0], np.rad2deg(theta[-1]), d[-1])
        diagnostics["log_hough"] = hv.Image(np.log(1 + h), bounds=bounds)
        # TODO: fix the left-edge vs. right-edge issue for theta bins
        diff_h_plot = hv.Curve((np.rad2deg(theta[:-1]), diff_h_std))
        if smooth:
            diff_h_plot *= hv.Curve(
                (np.rad2deg(theta[:-1]), diff_h_std_smoothed)
            ).options(color="cyan")
        diff_h_plot *= hv.VLine(np.rad2deg(angle)).options(color="red")
        diagnostics["diff_h_std"] = diff_h_plot
        diagnostics["angle"] = np.rad2deg(angle)
    return angle


def detect_rotation(img, window=np.deg2rad(10), diagnostics=None):
    angle1 = find_hough_angle(
        img, diagnostics=getitem_if_not_none(diagnostics, "hough_1")
    )
    angle2 = find_hough_angle(
        img,
        theta=np.linspace(angle1 - window, angle1 + window, 200),
        diagnostics=getitem_if_not_none(diagnostics, "hough_2"),
    )
    return np.pi / 2 - angle2


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


def detect_trench_set(img, theta, diagnostics=None):
    t0 = detect_trench_region(
        img, theta, diagnostics=getitem_if_not_none(diagnostics, "trench_region")
    )
    trench_anchors = detect_trench_anchors(
        img, t0, theta, diagnostics=getitem_if_not_none(diagnostics, "trench_anchors")
    )
    trench_points = detect_trench_ends(
        img,
        trench_anchors,
        theta,
        diagnostics=getitem_if_not_none(diagnostics, "trench_ends"),
    )
    return trench_points


def find_trench_threshold(img, bins=10, diagnostics=None):
    if diagnostics is not None:
        threshold_img_x = np.zeros(img.shape)
        threshold_img_y = np.zeros(img.shape)
    thresholds_x = []
    xs = np.linspace(0, img.shape[1], bins).astype(np.int_)
    for x0, x1 in zip(xs[:-1], xs[1:]):
        threshold_x = skimage.filters.threshold_otsu(img[:, x0:x1])
        if diagnostics is not None:
            threshold_img_x[:, x0:x1] = threshold_x
        thresholds_x.append(threshold_x)
    thresholds_y = []
    ys = np.linspace(0, img.shape[0], bins).astype(np.int_)
    for y0, y1 in zip(ys[:-1], ys[1:]):
        threshold_y = skimage.filters.threshold_otsu(img[y0:y1, :])
        if diagnostics is not None:
            threshold_img_y[y0:y1, :] = threshold_y
        thresholds_y.append(threshold_y)
    threshold = np.median(thresholds_x)
    if diagnostics is not None:
        diagnostics["threshold_img_x"] = RevImage(threshold_img_x)
        diagnostics["threshold_img_y"] = RevImage(threshold_img_y)
        diagnostics["threshold"] = threshold
    return threshold


def label_for_trenches(
    img,
    min_component_size=30,
    max_component_size=10**4,
    lowpass_radius=100,
    approximate_gaussian_filter=True,
    diagnostics=None,
):
    img = img.astype(np.float_)
    if diagnostics is not None:
        diagnostics["image"] = RevImage(img)
    if approximate_gaussian_filter:
        img_lowpass = gaussian_box_approximation(img, lowpass_radius)
    else:
        img_lowpass = skimage.filters.gaussian(img, lowpass_radius)
    img_highpass = img - img_lowpass
    if diagnostics is not None:
        diagnostics["lowpass_image"] = RevImage(img_lowpass)
        diagnostics["highpass_image"] = RevImage(img_highpass)
    threshold = find_trench_threshold(
        img_highpass,
        diagnostics=getitem_if_not_none(diagnostics, "find_trench_threshold"),
    )
    img_thresh = img_highpass > threshold
    if diagnostics is not None:
        diagnostics["thresholded_image"] = RevImage(img_thresh)
    # img_thresh = skimage.morphology.binary_erosion(img_thresh)
    # if diagnostics is not None:
    #     diagnostics['eroded_image'] = RevImage(img_thresh)
    components, num_components = skimage.morphology.label(img_thresh, return_num=True)
    if diagnostics is not None:
        diagnostics["components"] = RevImage(components)
        diagnostics["num_components"] = num_components
    cleaned_components = components.copy()
    skimage.morphology.remove_small_objects(
        cleaned_components, min_size=min_component_size, in_place=True
    )
    remove_large_objects(cleaned_components, max_component_size, in_place=True)
    cleaned_components, _, inverse_map = skimage.segmentation.relabel_sequential(
        cleaned_components
    )
    num_cleaned_components = len(inverse_map)
    if diagnostics is not None:
        diagnostics["cleaned_components"] = RevImage(cleaned_components)
        diagnostics["num_cleaned_components"] = num_cleaned_components
    normalized_img = normalize_componentwise(
        img_highpass,
        cleaned_components,
        label_index=np.arange(num_cleaned_components) + 1,
    )  # TODO: check arange
    if diagnostics is not None:
        diagnostics["normalized_image"] = RevImage(normalized_img)
    img_labels, label_index = label_binary_image(cleaned_components)
    if diagnostics is not None:
        diagnostics["labeled_image"] = RevImage(img_labels)
        diagnostics["label_index"] = tuple(label_index)
    return normalized_img, img_labels, label_index


def _get_trench_set(img, img_mask, theta=None, diagnostics=None):
    img_masked = np.where(
        skimage.morphology.binary_dilation(img_mask), img, np.percentile(img, 5)
    )
    # TODO: should only need to detect rotation once per image, not per trench set
    if theta is None:
        theta = detect_rotation(
            img_masked, diagnostics=getitem_if_not_none(diagnostics, "trench_rotation")
        )
    trench_points = detect_trench_set(img_masked, theta, diagnostics=diagnostics)
    return trench_points


def get_trenches_periodogram(img, setwise=True, diagnostics=None):
    normalized_img, img_labels, label_index = label_for_trenches(
        img, diagnostics=getitem_if_not_none(diagnostics, "labeling")
    )
    if setwise:
        angle = None
    else:
        angle = detect_rotation(
            img, diagnostics=getitem_if_not_none(diagnostics, "trench_rotation")
        )
    trenches = {}
    for label in label_index:
        trenches[label] = _get_trench_set(
            normalized_img,
            img_labels == label,
            angle=angle,
            diagnostics=getitem_if_not_none(diagnostics, "label_{}".format(label)),
        )
    return trenches
