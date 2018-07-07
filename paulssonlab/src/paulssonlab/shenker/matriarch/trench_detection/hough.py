import numpy as np
import scipy
import skimage
import skimage.morphology
import skimage.segmentation
from itertools import zip_longest
import holoviews as hv
import holoviews.operation.datashader as datashader

# TODO: fix imports
from holoborodko_diff import holo_diff
from .core import edge_point, coords_along, stack_jagged, stack_jagged_points
from .periodogram import label_for_trenches
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


def find_periodic_lines(
    img,
    theta=None,
    smooth=4,
    hough_func=hough_line_intensity,
    num_offset_points=200,
    diagnostics=None,
):
    if theta is None:
        theta = np.linspace(np.deg2rad(-45), np.deg2rad(45), 90)
    h, theta, d = hough_func(img, theta=theta)
    diff_h = np.diff(h.astype(np.int_), axis=1)  # TODO: is diff necessary??
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
        diagnostics["input"] = RevImage(img)
        diagnostics["angle_range"] = (np.rad2deg(theta[0]), np.rad2deg(theta[-1]))
        bounds = (np.rad2deg(theta[0]), d[0], np.rad2deg(theta[-1]), d[-1])
        diagnostics["log_hough"] = hv.Image(np.log(1 + h), bounds=bounds)
        # TODO: fix the left-edge vs. right-edge issue for theta bins
        theta_degrees = np.rad2deg(theta[:-1])
        diff_h_plot = hv.Curve((theta_degrees, diff_h_std))
        if smooth:
            diff_h_plot *= hv.Curve((theta_degrees, diff_h_std_smoothed)).options(
                color="cyan"
            )
        diff_h_plot *= hv.VLine(np.rad2deg(angle)).options(color="red")
        diagnostics["diff_h_std"] = diff_h_plot
        diagnostics["angle"] = np.rad2deg(angle)
    profile = h[:, theta_idx]
    freqs, spectrum = scipy.signal.periodogram(profile, scaling="spectrum")
    spectrum[:2] = 0  # TODO: ignore two highest frequencies
    pitch_idx = spectrum.argmax()
    pitch = 1 / freqs[pitch_idx]
    if diagnostics is not None:
        diagnostics["pitch"] = pitch
        diagnostics["spectrum"] = hv.Curve(spectrum) * hv.VLine(pitch_idx).options(
            color="red"
        )
    offsets = np.linspace(0, pitch, num_offset_points, endpoint=False)
    offset_idxs = (np.arange(0, len(profile), pitch) + offsets[:, np.newaxis]).astype(
        np.int_
    )
    offset_objective = profile[offset_idxs].sum(axis=1)
    if smooth:
        offset_objective_smoothed = scipy.ndimage.filters.gaussian_filter1d(
            offset_objective, smooth
        )
    else:
        offset_objective_smoothed = offset_objective
    offset_idx = offset_objective_smoothed.argmax()
    offset = offsets[offset_idx]
    if diagnostics is not None:
        diagnostics["offset"] = offset
        offset_plot = hv.Curve((offsets, offset_objective))
        if smooth:
            offset_plot *= hv.Curve((offsets, offset_objective_smoothed)).options(
                color="cyan"
            )
        offset_plot *= hv.VLine(offset).options(color="red")
        diagnostics["offsets"] = offset_plot
        highlighted_points = hv.Scatter(
            (d[offset_idxs[offset_idx]], profile[offset_idxs[offset_idx]])
        ).options(size=5, color="red")
        diagnostics["profile"] = hv.Curve((d, profile)) * highlighted_points
    return angle, pitch, offset


def find_trench_lines(img, window=np.deg2rad(10), diagnostics=None):
    angle1, pitch1, offset1 = find_periodic_lines(
        img, diagnostics=getitem_if_not_none(diagnostics, "hough_1")
    )
    angle2, pitch2, offset2 = find_periodic_lines(
        img,
        theta=np.linspace(angle1 - window, angle1 + window, 200),
        diagnostics=getitem_if_not_none(diagnostics, "hough_2"),
    )
    return angle2, pitch2, offset2


def trench_anchors(angle, pitch, offset, x_lim, y_lim):
    x_min, x_max = x_lim
    y_min, y_max = y_lim
    max_dist = int(np.ceil(np.sqrt(x_max**2 + y_max**2)))
    if angle < 0:
        effective_offset = (x_max * np.cos(angle) % pitch) - offset
    else:
        effective_offset = ((max_dist * 2) % pitch) + offset
    # TODO: the 2*pitch is just so offsets won't push rhos away from
    #       starting point
    # rhos = (np.arange(-2*pitch, x_max/np.cos(angle)+2*pitch, pitch) + effective_offset)
    # rhos = rhos[(rhos > 0) & (rhos < x_max/np.cos(angle))]
    abs_angle = np.abs(angle)
    delta = (y_max - x_max * np.tan(abs_angle)) * np.sin(abs_angle)
    rhos = np.arange(effective_offset % pitch, x_max / np.cos(angle) + delta, pitch)
    upper_right = np.array((x_max, 0))
    anchors = []
    for rho in rhos:
        anchor = rho * np.array((np.cos(angle), np.sin(angle)))
        if angle < 0:
            anchor = upper_right - anchor
        anchors.append(anchor)
    return anchors


def get_trench_profiles(img, anchors, angle, diagnostics=None):
    x_lim, y_lim = get_image_limits(img.shape)
    profile_points = []
    trench_profiles = []
    for anchor in anchors:
        x_end = edge_point(anchor, angle, x_lim, y_lim)
        xs, ys = coords_along(anchor, x_end)
        points = np.vstack((xs, ys)).T
        profile_points.append(points)
        trench_profiles.append(img[ys, xs])
    if diagnostics is not None:
        # TODO: make hv.Path??
        diagnostics["trench_profiles"] = hv.Overlay.from_values(
            [hv.Curve(tp) for tp in trench_profiles]
        )
    return trench_profiles, profile_points


def find_trench_ends(
    img, angle, pitch, offset, margin=15, threshold=0.8, smooth=100, diagnostics=None
):
    x_lim, y_lim = get_image_limits(img.shape)
    anchors = trench_anchors(angle, pitch, offset, x_lim, y_lim)
    top_profiles, top_profile_points = get_trench_profiles(
        img,
        anchors,
        3 / 2 * np.pi - angle,
        diagnostics=getitem_if_not_none(diagnostics, "top"),
    )
    bottom_profiles, bottom_profile_points = get_trench_profiles(
        img,
        anchors,
        np.pi / 2 - angle,
        diagnostics=getitem_if_not_none(diagnostics, "bottom"),
    )
    if diagnostics is not None:
        diagnostics["threshold"] = threshold
        lines_plot = hv.Path(
            [
                [top[-1], bottom[-1]]
                for top, bottom in zip(top_profile_points, bottom_profile_points)
            ]
        ).options(color="blue")
        top_line_plot = hv.Points([top[-1] for top in top_profile_points]).options(
            color="green"
        )
        bottom_line_plot = hv.Points(
            [bottom[-1] for bottom in bottom_profile_points]
        ).options(color="red")
        anchor_points_plot = hv.Points(anchors).options(size=3, color="cyan")
        diagnostics["image_with_lines"] = (
            RevImage(img)
            * lines_plot
            * top_line_plot
            * bottom_line_plot
            * anchor_points_plot
        )
    stacked_top_profile = np.nanpercentile(
        stack_jagged(top_profiles), threshold * 100, axis=0
    )
    stacked_bottom_profile = np.nanpercentile(
        stack_jagged(bottom_profiles), threshold * 100, axis=0
    )
    anchor_idx = len(stacked_top_profile)
    stacked_profile = np.concatenate(
        (stacked_top_profile[:0:-1], stacked_bottom_profile)
    )
    stacked_profile_diff = holo_diff(1, stacked_profile)
    top_end = max(stacked_profile_diff.argmax() - margin, 0)
    bottom_end = min(stacked_profile_diff.argmin() + margin, len(stacked_profile) - 1)
    if diagnostics is not None:
        diagnostics["margin"] = margin
        diagnostics["stacked_profile"] = (
            hv.Curve(stacked_profile)
            * hv.Curve(stacked_profile_diff).options(color="cyan")
            * hv.VLine(anchor_idx).options(color="gray")
            * hv.VLine(top_end).options(color="green")
            * hv.VLine(bottom_end).options(color="red")
        )
    stacked_top_points = stack_jagged_points(top_profile_points)
    stacked_bottom_points = stack_jagged_points(bottom_profile_points)
    stacked_points = np.concatenate((stacked_top_points[:0:-1], stacked_bottom_points))
    top_endpoints = stacked_points[top_end]
    bottom_endpoints = stacked_points[bottom_end]
    print("a", anchors[-1], len(top_profile_points))
    print("at", top_profile_points[-1])
    print("ab", bottom_profile_points[-1])
    print("s", stacked_points.shape, top_end, bottom_end)
    print("stacked_points", stacked_points[top_end - 2 : top_end + 2])
    print("hhh", stacked_points[top_end - 100 : bottom_end + 100, -1])
    print("top", top_endpoints)
    # discard trenches where top endpoint is the same as the bottom endpoint
    mask = ~np.apply_along_axis(np.all, 1, np.equal(top_endpoints, bottom_endpoints))
    top_endpoints = top_endpoints[mask]
    bottom_endpoints = bottom_endpoints[mask]
    if diagnostics is not None:
        trench_plot = hv.Path(
            [
                [top_endpoint, bottom_endpoint]
                for top_endpoint, bottom_endpoint in zip(
                    top_endpoints, bottom_endpoints
                )
            ]
        ).options(color="white")
        top_points_plot = hv.Points(top_endpoints).options(size=3, color="green")
        bottom_points_plot = hv.Points(bottom_endpoints).options(size=3, color="red")
        diagnostics["image_with_trenches"] = (
            RevImage(img) * trench_plot * top_points_plot * bottom_points_plot
        )
    return top_endpoints, bottom_endpoints


def find_trenches(img, setwise=True, diagnostics=None):
    img_normalized, img_labels, label_index = label_for_trenches(
        img, diagnostics=getitem_if_not_none(diagnostics, "labeling")
    )
    trenches = {}
    if not setwise:
        img_masked = np.where(
            skimage.morphology.binary_dilation(img_labels != 0),
            img_normalized,
            np.percentile(img_normalized, 5),
        )
        angle, pitch, offset = find_trench_lines(
            img_masked,
            diagnostics=getitem_if_not_none(diagnostics, "find_trench_lines"),
        )
    for label in label_index:
        label_diagnostics = getitem_if_not_none(diagnostics, "label_{}".format(label))
        img_masked = np.where(
            skimage.morphology.binary_dilation(img_labels == label),
            img_normalized,
            np.percentile(img_normalized, 5),
        )
        if setwise:
            angle, pitch, offset = find_trench_lines(
                img_masked,
                diagnostics=getitem_if_not_none(label_diagnostics, "find_trench_lines"),
            )
        trenches[label] = find_trench_ends(
            img_masked,
            angle,
            pitch,
            offset,
            diagnostics=getitem_if_not_none(label_diagnostics, "find_trench_ends"),
        )
    return trenches
