import numpy as np
import pandas as pd
import holoviews as hv
from holoborodko_diff import holo_diff
from geometry import get_image_limits
from ui import RevImage
from workflow import points_dataframe
from .hough import trench_anchors
from .geometry import edge_point, coords_along


def get_trench_line_profiles(
    img, angle, anchor_rho, rho_min, rho_max, diagnostics=None
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
    if diagnostics is not None:
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
    profiles = np.array(padded_profiles)
    stacked_points = np.array(padded_line_points).swapaxes(0, 1)
    return profiles, stacked_points, anchor_idx


def find_trench_ends(
    img,
    angle,
    anchor_rho,
    rho_min,
    rho_max,
    margin=15,
    profile_quantile=0.95,
    diagnostics=None,
):
    profiles, stacked_points, anchor_idx = get_trench_line_profiles(
        img, angle, anchor_rho, rho_min, rho_max, diagnostics=diagnostics
    )
    stacked_profile = np.nanpercentile(profiles, profile_quantile * 100, axis=0)
    stacked_profile_diff = holo_diff(1, stacked_profile)
    # using np.nanargmax/min because we might have an all-nan axis
    top_end = max(np.nanargmax(stacked_profile_diff) - margin, 0)
    bottom_end = min(
        np.nanargmin(stacked_profile_diff) + margin, len(stacked_profile) - 1
    )
    if diagnostics is not None:
        diagnostics["profile_quantile"] = profile_quantile
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
