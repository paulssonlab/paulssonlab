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
