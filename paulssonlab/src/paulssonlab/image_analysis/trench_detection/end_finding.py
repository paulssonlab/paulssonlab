from collections import defaultdict

import holoviews as hv
import numpy as np
import pandas as pd
import scipy
import skimage.measure
import skimage.morphology

from paulssonlab.image_analysis.trench_detection.profile import angled_profiles
from paulssonlab.image_analysis.ui import RevImage
from paulssonlab.util.numeric import silent_nanquantile


def find_trench_ends(
    img,
    angle,
    rhos,
    profile_quantile=0.95,
    margin=5,
    min_length=50,
    smooth=10,
    threshold=0.2,
    threshold_quantile=0.99,
    diagnostics=None,
):
    profiles, profile_points = angled_profiles(
        img, angle, rhos, diagnostics=diagnostics
    )
    profiles = profiles.astype(np.float32)
    # treat <=0 values as background
    # (useful for uint16 images where we can't use NaN)
    profiles[profiles <= 0] = np.nan
    reduced_profile = silent_nanquantile(profiles, profile_quantile, axis=0)
    if smooth:
        reduced_profile_smooth = scipy.ndimage.filters.gaussian_filter1d(
            reduced_profile, smooth
        )
    else:
        reduced_profile_smooth = reduced_profile
    threshold_value = (
        silent_nanquantile(reduced_profile_smooth, threshold_quantile) * threshold
    )
    profile_mask = reduced_profile_smooth > threshold_value
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
        diagnostics["margin"] = margin
        start_lines = hv.Overlay([hv.VLine(x[0]) for x in endpoints]).options(
            "VLine", color="gray"
        ) * hv.Overlay([hv.VLine(max(x[0] - margin, 0)) for x in endpoints]).options(
            "VLine", color="green"
        )
        stop_lines = hv.Overlay([hv.VLine(x[1]) for x in endpoints]).options(
            "VLine", color="gray"
        ) * hv.Overlay(
            [hv.VLine(min(x[1] + margin, len(profile_points) - 1)) for x in endpoints]
        ).options(
            "VLine", color="red"
        )
        diagnostics["reduced_profile"] = (
            hv.Curve(reduced_profile)
            * hv.HLine(threshold_value).options(color="gray")
            * start_lines
            * stop_lines
        )
    trench_dfs = {}
    if diagnostics is not None:
        trench_lines = []
    for trench_set_idx, (top_end, bottom_end) in enumerate(endpoints):
        top_endpoints = profile_points[max(top_end - margin, 0)]
        bottom_endpoints = profile_points[
            min(bottom_end + margin, len(profile_points) - 1)
        ]
        # discard trenches where top endpoint is the same as the bottom endpoint
        # this also throws out out-of-range rhos which have their endpoints set to (nan, nan)
        mask = ~np.apply_along_axis(
            np.all, 1, np.isclose(top_endpoints, bottom_endpoints, equal_nan=True)
        )
        top_endpoints = top_endpoints[mask]
        bottom_endpoints = bottom_endpoints[mask]
        trench_idxs = np.arange(len(mask))[mask]
        trench_dfs[trench_set_idx] = pd.DataFrame(
            {
                "top_x": top_endpoints[:, 0],
                "top_y": top_endpoints[:, 1],
                "bottom_x": bottom_endpoints[:, 0],
                "bottom_y": bottom_endpoints[:, 1],
            },
            index=trench_idxs,
        ).rename_axis(index="trench_line")
        if diagnostics is not None:
            top_endpoints_shifted = top_endpoints + 0.5
            bottom_endpoints_shifted = bottom_endpoints + 0.5
            trench_lines.extend(
                [
                    [top_endpoint, bottom_endpoint]
                    for top_endpoint, bottom_endpoint in zip(
                        top_endpoints_shifted, bottom_endpoints_shifted
                    )
                ]
            )
    if diagnostics is not None:
        trench_plot = hv.Path(trench_lines).options(color="white")
        top_points_plot = hv.Points([line[0] for line in trench_lines]).options(
            size=3, color="green"
        )
        bottom_points_plot = hv.Points([line[1] for line in trench_lines]).options(
            size=3, color="red"
        )
        diagnostics["image_with_trenches"] = (
            RevImage(img) * trench_plot * top_points_plot * bottom_points_plot
        )
    df = pd.concat(trench_dfs, names=["trench_set"])
    return df
