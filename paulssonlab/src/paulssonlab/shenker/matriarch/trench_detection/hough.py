import numpy as np
import scipy
import skimage
import skimage.morphology
import skimage.segmentation
from itertools import zip_longest
import holoviews as hv
import holoviews.operation.datashader as datashader

# TODO: fix imports
from .core import detect_trench_ends
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
            (offset_idxs[offset_idx], profile[offset_idxs[offset_idx]])
        ).options(size=5, color="red")
        diagnostics["profile"] = hv.Curve(profile) * highlighted_points
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
    return np.pi / 2 - angle2, pitch2, offset2


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
        # trenches[label] = find_trench_ends(img_masked, diagnostics=getitem_if_not_none(label_diagnostics, 'find_trench_ends'))
    return trenches
