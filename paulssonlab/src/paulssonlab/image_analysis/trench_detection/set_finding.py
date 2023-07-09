from collections import defaultdict

import holoviews as hv
import numpy as np
import scipy
import skimage
import sklearn
import sklearn.cluster
from matplotlib.path import Path
from scipy.spatial import ConvexHull
from sklearn.preprocessing import StandardScaler

from paulssonlab.image_analysis.image import (
    gaussian_box_approximation,
    normalize_componentwise,
    remove_large_objects,
)
from paulssonlab.image_analysis.misc.holoborodko_diff import holo_diff
from paulssonlab.image_analysis.trench_detection.profile import get_trench_line_profiles
from paulssonlab.image_analysis.ui import RevImage
from paulssonlab.image_analysis.util import getitem_if_not_none
from paulssonlab.util.numeric import silent_nanquantile


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


def drop_rare_labels(labels):
    counter = Counter(labels)
    total = sum(counter)
    good_labels = []
    for label, count in counter.iteritems():
        if count / total > 0.01:
            good_labels.append(label)
    return good_labels


# TODO: remove
def find_trench_sets_by_clustering(
    img, img_mask, angle, anchor_rho, rho_min, rho_max, diagnostics=None
):
    X, fit = cluster_binary_image(img_mask)
    img_labels = np.zeros_like(img_mask, dtype=np.int8)  # TODO: fixing dtype
    for i in range(len(fit.labels_)):
        img_labels[X[i, 0], X[i, 1]] = fit.labels_[i] + 1
    label_index = np.sort(np.unique(fit.labels_)) + 1
    if diagnostics is not None:
        diagnostics["labeled_image"] = RevImage(img_labels)
        # diagnostics['label_index'] = tuple(label_index) # TODO: arrow/parquet nested column
        diagnostics["label_min"] = label_index[0]
        diagnostics["label_max"] = label_index[-1]
    return img_labels, label_index


def find_trench_sets_by_cutting(
    img,
    img_mask,
    angle,
    rhos,
    profile_quantile=0.95,
    min_length=50,
    smooth=10,
    threshold=0.05,
    diagnostics=None,
):
    profiles, stacked_points = get_trench_line_profiles(
        img, angle, rhos, diagnostics=diagnostics
    )
    stacked_profile = silent_nanquantile(profiles, profile_quantile, axis=0)
    if smooth:
        stacked_profile_smooth = scipy.ndimage.filters.gaussian_filter1d(
            stacked_profile, smooth
        )
    else:
        stacked_profile_smooth = stacked_profile
    threshold_value = silent_nanquantile(stacked_profile_smooth, threshold)
    profile_mask = stacked_profile_smooth > threshold_value
    if min_length:
        skimage.morphology.remove_small_objects(
            profile_mask, min_size=min_length, out=profile_mask
        )
    profile_labels = skimage.measure.label(profile_mask)
    endpoints = []
    for label in range(1, profile_labels.max() + 1):
        nonzero = np.nonzero(profile_labels == label)[0]
        endpoints.append((nonzero[0], nonzero[-1]))
    if diagnostics is not None:
        diagnostics["profile_quantile"] = profile_quantile
        diagnostics["threshold"] = threshold
        diagnostics["threshold_value"] = threshold_value
        start_lines = hv.Overlay([hv.VLine(x[0]) for x in endpoints]).options(
            "VLine", color="green"
        )
        stop_lines = hv.Overlay([hv.VLine(x[1]) for x in endpoints]).options(
            "VLine", color="red"
        )
        diagnostics["stacked_profile"] = (
            hv.Curve(stacked_profile)
            * hv.HLine(threshold_value).options(color="gray")
            * start_lines
            * stop_lines
        )
    label_index = tuple(range(1, len(endpoints) + 1))
    img_labels = np.zeros_like(img_mask, dtype=np.int8)
    width, height = img.shape
    y, x = np.nonzero(img_mask)
    coors = np.hstack((x.reshape(-1, 1), y.reshape(-1, 1)))
    for label, (top_end, bottom_end) in zip(label_index, endpoints):
        top_endpoints = stacked_points[top_end]
        bottom_endpoints = stacked_points[bottom_end]
        # discard trenches where top endpoint is the same as the bottom endpoint
        mask = ~np.apply_along_axis(
            np.all, 1, np.equal(top_endpoints, bottom_endpoints)
        )
        top_endpoints = top_endpoints[mask]
        bottom_endpoints = bottom_endpoints[mask]
        polygon = np.vstack((top_endpoints, bottom_endpoints))
        hull = ConvexHull(polygon)
        poly_path = Path(polygon[hull.vertices])
        # FROM: https://stackoverflow.com/questions/3654289/scipy-create-2d-polygon-mask
        mask = poly_path.contains_points(coors)
        img_labels[y[mask], x[mask]] = label
    if diagnostics is not None:
        diagnostics["labeled_image"] = RevImage(img_labels)
        diagnostics["label_min"] = label_index[0]
        diagnostics["label_max"] = label_index[-1]
    return img_labels, label_index


# TODO: WIP, remove?
def find_trench_sets_by_diff_cutting(
    img,
    img_mask,
    angle,
    anchor_rho,
    rho_min,
    rho_max,
    profile_quantile=0.95,
    min_length=50,
    smooth=10,
    prominence=0.5,
    diagnostics=None,
):
    profiles, _, _ = get_trench_line_profiles(
        img, angle, anchor_rho, rho_min, rho_max, diagnostics=diagnostics
    )
    stacked_profile = silent_nanquantile(profiles, profile_quantile, axis=0)
    if smooth:
        stacked_profile_smooth = scipy.ndimage.filters.gaussian_filter1d(
            stacked_profile, smooth
        )
    else:
        stacked_profile_smooth = stacked_profile
    stacked_profile_diff = holo_diff(1, stacked_profile_smooth)
    stacked_profile_diff[np.isnan(stacked_profile_diff)] = 0
    abs_prominence = np.abs(stacked_profile_diff).max() * prominence
    stacked_profile_diff_pos = np.clip(stacked_profile_diff, 0, None)
    stacked_profile_diff_neg = np.clip(-stacked_profile_diff, 0, None)
    start_idxs = scipy.signal.find_peaks(
        stacked_profile_diff_pos,
        distance=min_length,
        prominence=abs_prominence,
    )
    stop_idxs = scipy.signal.find_peaks(
        stacked_profile_diff_neg,
        distance=min_length,
        prominence=abs_prominence,
    )
    if diagnostics is not None:
        diagnostics["profile_quantile"] = profile_quantile
        start_lines = hv.Overlay([hv.VLine(x) for x in start_idxs]).options(
            "VLine", color="green"
        )
        stop_lines = hv.Overlay([hv.VLine(x) for x in stop_idxs]).options(
            "VLine", color="red"
        )
        diagnostics["stacked_profile"] = (
            hv.Curve(stacked_profile)
            * hv.Curve(stacked_profile_diff).options(color="cyan")
            * start_lines
            * stop_lines
        )
    #             hv.VLine(start_idxs).options(color='green') * \
    #             hv.VLine(stop_idxs).options(color='red')
    #     if diagnostics is not None:
    #         diagnostics['labeled_image'] = RevImage(img_labels)
    #         diagnostics['label_index'] = tuple(label_index)
    raise NotImplementedError
    return find_trench_sets_by_clustering(
        img, img_mask, angle, anchor_rho, rho_min, rho_max, diagnostics=diagnostics
    )
    return img_labels, label_index


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


def binarize_trench_image(
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
        cleaned_components, min_size=min_component_size, out=cleaned_components
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
        label_index=np.arange(num_cleaned_components) + 1,  # TODO: is this right?
    )  # TODO: check arange
    if diagnostics is not None:
        diagnostics["normalized_image"] = RevImage(normalized_img)
    return normalized_img, cleaned_components
