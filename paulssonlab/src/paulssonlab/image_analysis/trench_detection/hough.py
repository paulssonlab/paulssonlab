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
from ..misc.holoborodko_diff import holo_diff
import peakutils
from .geometry import edge_point, coords_along
from .peaks import find_periodic_peaks
from ..util import getitem_if_not_none
from ..ui import RevImage
from ..image import (
    hough_line_intensity,
    remove_large_objects,
    normalize_componentwise,
    gaussian_box_approximation,
)
from ..geometry import get_image_limits
from ..workflow import points_dataframe
from .. import common


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


def find_peaks(profile, threshold=0.2, min_dist=5, diagnostics=None):
    idxs = peakutils.indexes(profile, thres=threshold, min_dist=min_dist)
    return idxs, None


def find_periodic_lines(
    img,
    theta=None,
    smooth=4,
    threshold_quantile=0.75,
    upscale=None,
    hough_func=hough_line_intensity,
    peak_func=find_periodic_peaks,
    diagnostics=None,
):
    if theta is None:
        theta = np.linspace(np.deg2rad(-45), np.deg2rad(45), 90)
    if threshold_quantile is not None:
        img_thresh = img * (img > np.percentile(img, threshold_quantile * 100))
    else:
        img_thresh = img
    h, theta, rho = hough_func(img_thresh, theta=theta)
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
        diagnostics["thresholded_image"] = RevImage(img_thresh)
        # diagnostics['angle_range'] = (np.rad2deg(theta[0]), np.rad2deg(theta[-1])) # TODO: arrow/parquet nested column
        diagnostics["angle_range_min"] = np.rad2deg(theta[0])
        diagnostics["angle_range_max"] = np.rad2deg(theta[-1])
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
    # TODO: can probably reduce sampling here in one or both hough calls
    angle1, *_ = find_periodic_lines(
        img, diagnostics=getitem_if_not_none(diagnostics, "hough_1")
    )
    res2 = find_periodic_lines(
        img,
        theta=np.linspace(angle1 - window, angle1 + window, 200),
        diagnostics=getitem_if_not_none(diagnostics, "hough_2"),
    )
    return res2
