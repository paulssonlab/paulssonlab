import numpy as np
import pandas as pd
import scipy
import skimage
import skimage.morphology
import skimage.segmentation
from itertools import zip_longest
import holoviews as hv
import holoviews.operation.datashader as datashader

# TODO: fix imports
from holoborodko_diff import holo_diff
import peakutils
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
from workflow import points_dataframe
import common


def find_peaks(profile, threshold=0.2, min_dist=5, diagnostics=None):
    idxs = peakutils.indexes(profile, thres=threshold, min_dist=min_dist)
    return idxs, None


def find_periodic_peaks(
    profile,
    refine=True,
    nfft=2**14,
    smooth_offset=4,
    num_offset_points=200,
    diagnostics=None,
):
    freqs, spectrum = scipy.signal.periodogram(
        profile, window="hann", nfft=nfft, scaling="spectrum"
    )
    pitch_idx = spectrum.argmax()
    pitch = 1 / freqs[pitch_idx]
    if diagnostics is not None:
        diagnostics["pitch"] = pitch
        spectrum_plot = hv.Curve(spectrum)
        spectrum_plot *= hv.VLine(pitch_idx).options(color="red")
        diagnostics["spectrum"] = spectrum_plot
    offsets = np.linspace(0, pitch, num_offset_points, endpoint=False)
    # always leave a period of length `pitch` so we can add offsets
    # even when we could've fit another period
    # TODO: this sometimes results in one trench undetected!
    offset_idxs = (
        np.arange(0, len(profile) - pitch, pitch) + offsets[:, np.newaxis]
    ).astype(np.int_)
    offset_objective = profile[offset_idxs].sum(axis=1)
    if smooth_offset:
        offset_objective_smoothed = scipy.ndimage.filters.gaussian_filter1d(
            offset_objective, smooth_offset
        )
    else:
        offset_objective_smoothed = offset_objective
    offset_idx = offset_objective_smoothed.argmax()
    offset = offsets[offset_idx]
    idxs = offset_idxs[offset_idx]
    if diagnostics is not None:
        diagnostics["offset"] = offset
        offset_plot = hv.Curve((offsets, offset_objective))
        if smooth_offset:
            offset_plot *= hv.Curve((offsets, offset_objective_smoothed)).options(
                color="cyan"
            )
        offset_plot *= hv.VLine(offset).options(color="red")
        diagnostics["offsets"] = offset_plot
    if not refine:
        refined_idxs = idxs
        info = None
    else:
        idx_start = np.clip(idxs - refine, 0, len(profile))
        idx_end = np.clip(idxs + refine, 0, len(profile))
        shifts = np.array(
            [
                profile[idx_start[i] : idx_end[i]].argmax() - (idxs[i] - idx_start[i])
                for i in range(len(idxs))
            ]
        )
        shifts = np.where(profile[idxs] != profile[idxs + shifts], shifts, 0)
        refined_idxs = idxs + shifts
        info = pd.DataFrame({"hough_unshifted_value": profile[idxs], "shift": shifts})
        if diagnostics is not None:
            periodic_points = hv.Scatter((idxs, profile[idxs])).options(
                size=5, color="red"
            )
            refined_points = hv.Scatter((refined_idxs, profile[refined_idxs])).options(
                size=3, color="cyan"
            )
            diagnostics["refined_points"] = (
                hv.Curve(profile) * periodic_points * refined_points
            )
    return refined_idxs, info


def find_periodic_lines(
    img,
    theta=None,
    smooth=4,
    upscale=None,
    hough_func=hough_line_intensity,
    peak_func=find_periodic_peaks,
    diagnostics=None,
):
    if theta is None:
        theta = np.linspace(np.deg2rad(-45), np.deg2rad(45), 90)
    h, theta, rho = hough_func(img, theta=theta)
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
        bounds = (np.rad2deg(theta[0]), rho[0], np.rad2deg(theta[-1]), rho[-1])
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
    x_lim, y_lim = get_image_limits(img.shape)
    x_min, x_max = x_lim
    y_min, y_max = y_lim
    max_dist = int(np.ceil(np.sqrt(x_max**2 + y_max**2)))
    if angle < 0:
        rho_min = y_max * np.sin(angle)
        rho_max = x_max * np.cos(angle)
    else:
        rho_min = 0
        rho_max = max_dist * np.cos(angle - np.pi / 4)
    idx_min = np.where(rho_min <= rho)[0][0]
    idx_max = np.where(rho_max > rho)[0][-1]
    assert rho_min <= rho[idx_min] >= rho_min
    assert rho_min <= rho[idx_max] <= rho_max
    trimmed_profile = profile[idx_min : idx_max + 1]
    trimmed_rho = rho[idx_min : idx_max + 1]
    trimmed_profile_plot = hv.Curve((trimmed_rho, trimmed_profile))
    if upscale:
        interpolated_func = scipy.interpolate.interp1d(
            trimmed_rho, trimmed_profile, kind="cubic", copy=False, assume_sorted=True
        )
        interpolated_rho = np.linspace(
            trimmed_rho[0], trimmed_rho[-1], len(trimmed_rho) * upscale
        )
        interpolated_profile = interpolated_func(interpolated_rho)
        trimmed_profile_plot *= hv.Curve(
            (interpolated_rho, interpolated_profile)
        ).options(color="cyan")
    else:
        interpolated_profile = trimmed_profile
    if diagnostics is not None:
        diagnostics["trimmed_profile"] = trimmed_profile_plot
    anchor_idxs, anchor_info = peak_func(
        interpolated_profile, diagnostics=getitem_if_not_none(diagnostics, "peak_func")
    )
    if upscale:
        anchor_idxs = (anchor_idxs / upscale).astype(np.int_)
    anchor_values = trimmed_profile[anchor_idxs]  # TODO: interpolation
    anchor_idxs += idx_min
    anchor_rho = rho[anchor_idxs]
    # TODO interpolate rho
    if diagnostics is not None:
        rho_points = hv.Scatter((anchor_rho, anchor_values)).options(
            size=5, color="cyan"
        )
        profile_plot = hv.Curve((rho, profile)) * rho_points
        profile_plot *= hv.VLine(rho_min).options(color="red")
        profile_plot *= hv.VLine(rho_max).options(color="red")
        diagnostics["profile"] = profile_plot
    info = pd.DataFrame({"hough_value": anchor_values})
    if anchor_info is not None:
        info = info.join(anchor_info)
    return angle, anchor_rho, rho_min, rho_max, info


def find_trench_lines(img, window=np.deg2rad(10), diagnostics=None):
    angle1, *_ = find_periodic_lines(
        img, diagnostics=getitem_if_not_none(diagnostics, "hough_1")
    )
    res2 = find_periodic_lines(
        img,
        theta=np.linspace(angle1 - window, angle1 + window, 200),
        diagnostics=getitem_if_not_none(diagnostics, "hough_2"),
    )
    return res2


def trench_anchors(angle, anchor_rho, rho_min, rho_max, x_lim, y_lim):
    x_min, x_max = x_lim
    y_min, y_max = y_lim
    if angle < 0:
        anchor_rho = rho_max - anchor_rho
    anchors = (
        anchor_rho[:, np.newaxis]
        * np.array((np.cos(angle), np.sin(angle)))[np.newaxis, :]
    )
    if angle < 0:
        upper_right = np.array((x_max, 0))
        anchors = upper_right - anchors
    return anchors


def find_trench_ends(
    img,
    angle,
    anchor_rho,
    rho_min,
    rho_max,
    margin=15,
    threshold=0.95,
    smooth=100,
    diagnostics=None,
):
    x_lim, y_lim = get_image_limits(img.shape)
    anchors = trench_anchors(angle, anchor_rho, rho_min, rho_max, x_lim, y_lim)
    profiles = []
    line_points = []
    offsets = []
    for anchor in anchors:
        top_anchor = edge_point(anchor, 3 / 2 * np.pi - angle, x_lim, y_lim)
        bottom_anchor = edge_point(anchor, np.pi / 2 - angle, x_lim, y_lim)
        line_length = np.linalg.norm(top_anchor - bottom_anchor)
        top_length = np.linalg.norm(top_anchor - anchor)
        bottom_length = np.linalg.norm(bottom_anchor - anchor)
        xs, ys = coords_along(top_anchor, bottom_anchor)
        profile = img[ys, xs]
        points = np.vstack((xs, ys)).T
        # TODO: precision??
        if line_length >= max(top_length, bottom_length):
            # line contains anchor
            offset = -int(np.ceil(top_length))
        else:
            # line is strictly on one side of anchor
            if top_length <= bottom_length:
                # line lies below anchor
                offset = int(np.ceil(top_length))
            else:
                # line lies above anchor
                offset = int(np.ceil(bottom_length))
        profiles.append(profile)
        line_points.append(points)
        offsets.append(offset)
    min_offset = min(offsets)
    anchor_idx = -min_offset
    max_stacked_length = max(
        [len(profile) - offset for offset, profile in zip(offsets, profiles)]
    )
    padded_profiles = []
    padded_line_points = []
    for i, (profile, points, offset) in enumerate(zip(profiles, line_points, offsets)):
        left_padding = offset - min_offset
        right_padding = max_stacked_length - left_padding - len(profile)
        padded_profile = np.pad(
            profile, (left_padding, right_padding), "constant", constant_values=np.nan
        )
        padded_points = np.pad(points, [(left_padding, right_padding), (0, 0)], "edge")
        padded_profiles.append(padded_profile)
        padded_line_points.append(padded_points)
    if diagnostics is not None:
        # TODO: make hv.Path??
        diagnostics["profiles"] = hv.Overlay.from_values(
            [hv.Curve(tp) for tp in padded_profiles]
        )
    # stacked_profile = padded_profiles[25]
    stacked_profile = np.nanpercentile(
        np.array(padded_profiles), threshold * 100, axis=0
    )
    stacked_points = np.array(padded_line_points).swapaxes(0, 1)
    if diagnostics is not None:
        diagnostics["threshold"] = threshold
        lines_plot = hv.Path(
            [[points[0], points[-1]] for points in line_points]
        ).options(color="blue")
        top_line_plot = hv.Points([points[0] for points in line_points]).options(
            color="green"
        )
        bottom_line_plot = hv.Points([points[-1] for points in line_points]).options(
            color="red"
        )
        anchor_points_plot = hv.Points(anchors).options(size=3, color="cyan")
        diagnostics["image_with_lines"] = (
            RevImage(img)
            * lines_plot
            * top_line_plot
            * bottom_line_plot
            * anchor_points_plot
        )
    stacked_profile_diff = holo_diff(1, stacked_profile)
    # using np.nanargmax/min because we might have an all-nan axis
    top_end = max(np.nanargmax(stacked_profile_diff) - margin, 0)
    bottom_end = min(
        np.nanargmin(stacked_profile_diff) + margin, len(stacked_profile) - 1
    )
    if diagnostics is not None:
        diagnostics["margin"] = margin
        diagnostics["stacked_profile"] = (
            hv.Curve(stacked_profile)
            * hv.Curve(stacked_profile_diff).options(color="cyan")
            * hv.VLine(anchor_idx).options(color="gray")
            * hv.VLine(top_end).options(color="green")
            * hv.VLine(bottom_end).options(color="red")
        )
    top_endpoints = stacked_points[top_end]
    bottom_endpoints = stacked_points[bottom_end]
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
    trench_index = np.arange(len(mask))[mask]
    df_columns = {
        "top": points_dataframe(top_endpoints, index=trench_index),
        "bottom": points_dataframe(bottom_endpoints, index=trench_index),
    }
    df = pd.concat(df_columns, axis=1)
    return df


def find_trenches(img, reindex=True, diagnostics=None):
    img_normalized, img_labels, label_index = label_for_trenches(
        img, diagnostics=getitem_if_not_none(diagnostics, "labeling")
    )
    trench_sets = {}
    for label in label_index:
        label_diagnostics = getitem_if_not_none(diagnostics, "label_{}".format(label))
        img_masked = np.where(
            skimage.morphology.binary_dilation(img_labels == label),
            img_normalized,
            np.percentile(img_normalized, 5),
        )
        angle, anchor_rho, rho_min, rho_max, anchor_info = find_trench_lines(
            img_masked,
            diagnostics=getitem_if_not_none(label_diagnostics, "find_trench_lines"),
        )
        trench_sets[label] = find_trench_ends(
            img_masked,
            angle,
            anchor_rho,
            rho_min,
            rho_max,
            diagnostics=getitem_if_not_none(label_diagnostics, "find_trench_ends"),
        )
        if anchor_info is not None:
            anchor_info.columns = [("info", col) for col in anchor_info.columns]
            trench_sets[label] = trench_sets[label].join(anchor_info, how="left")
            if reindex:
                trench_sets[label].reset_index(drop=True, inplace=True)
    trenches_df = pd.concat(trench_sets)
    trenches_df.index.names = ["trench_set", "trench"]
    return trenches_df
