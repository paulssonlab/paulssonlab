import holoviews as hv
import numpy as np
import pandas as pd
from skimage.measure.profile import _line_profile_coordinates, profile_line
from skspatial.objects import Line, LineSegment

from paulssonlab.image_analysis.geometry import get_image_limits
from paulssonlab.image_analysis.image import hough_bounds
from paulssonlab.image_analysis.misc.holoborodko_diff import holo_diff
from paulssonlab.image_analysis.trench_detection.geometry import (
    edge_point,
    intersect_line_with_segment,
    trench_anchors,
)
from paulssonlab.image_analysis.ui import RevImage
from paulssonlab.image_analysis.workflow import points_dataframe
from paulssonlab.util.numeric import silent_nanquantile


def angled_line_profile_endpoints(angle, rho):
    # the top/bottom labels correspond to the coördinate system
    # where the y-axis is inverted (e.g., what RevImage assumes)
    top_line = LineSegment((x_lim[0], y_lim[0]), (x_lim[1], y_lim[0]))
    bottom_line = LineSegment((x_lim[0], y_lim[1]), (x_lim[1], y_lim[1]))
    left_line = LineSegment((x_lim[0], y_lim[0]), (x_lim[0], y_lim[1]))
    right_line = LineSegment((x_lim[1], y_lim[0]), (x_lim[1], y_lim[1]))
    anchor = (rho * np.cos(angle), 0)
    points = []
    for bounding_line in (top_line, bottom_line, left_line, right_line):
        intersection = intersect_line_with_segment(
            angled_line(anchor, angle), bounding_line
        )
        if intersection is not None:
            points.append(intersection)
    points = list(sorted(points, key=lambda x: x.distance_point(anchor)))
    if not points:
        return None, None, None
    assert len(points) == 2
    top, bottom = points
    if angle <= 0:
        corner = (x_lim[0], y_lim[0])
    else:
        corner = (x_lim[1], y_lim[0])
    offset = angled_line(anchor, angle + np.pi / 2).distance_point(top)
    return top, bottom, offset


def angled_profiles(img, angle, rhos, diagnostics=None):
    x_lim, y_lim = get_image_limits(img.shape)
    profiles = []
    line_points = []
    offsets = []
    for rho in rhos:
        top, bottom, offset = angled_line_profile_endpoints(angle, rho)
        assert top is not None  # if rho bounds are correct
        # need to give coordinates in (y, x) order
        profile = profile_line(
            img, top[::-1], bottom[::-1], mode="constant", cval=np.nan
        )
        points = _line_profile_coordinates(top[::-1], bottom[::-1]).swapaxes(0, 1)[
            :, ::-1, 0
        ]
        profiles.append(profile)
        line_points.append(points)
        offsets.append(offset)
    min_offset = min(offsets)
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
        diagnostics["profiles"] = hv.Overlay([hv.Curve(tp) for tp in padded_profiles])
    if diagnostics is not None:
        lines_plot = hv.Path(
            [[points[0] + 0.5, points[-1] + 0.5] for points in line_points]
        ).options(color="blue")
        top_line_plot = hv.Points([points[0] + 0.5 for points in line_points]).options(
            color="green"
        )
        bottom_line_plot = hv.Points(
            [points[-1] + 0.5 for points in line_points]
        ).options(color="red")
        anchor_points_plot = hv.Points(anchors + 0.5).options(size=3, color="cyan")
        diagnostics["image_with_lines"] = (
            RevImage(img)
            * lines_plot
            * top_line_plot
            * bottom_line_plot
            * anchor_points_plot
        )
    profiles = np.array(padded_profiles)
    stacked_points = np.array(padded_line_points).swapaxes(0, 1)
    return profiles, stacked_points


def get_trench_line_profiles(img, angle, rhos, diagnostics=None):
    x_lim, y_lim = get_image_limits(img.shape)
    rho_min, rho_max = hough_bounds(img.shape, angle)
    anchors = trench_anchors(angle, rhos, rho_min, rho_max, x_lim, y_lim)
    profiles = []
    line_points = []
    offsets = []
    for anchor in anchors:
        top_anchor = edge_point(anchor, 3 / 2 * np.pi - angle, x_lim, y_lim)
        bottom_anchor = edge_point(anchor, np.pi / 2 - angle, x_lim, y_lim)
        line_length = np.linalg.norm(top_anchor - bottom_anchor)
        top_length = np.linalg.norm(top_anchor - anchor)
        bottom_length = np.linalg.norm(bottom_anchor - anchor)
        # need to give coordinates in (y, x) order
        profile = profile_line(
            img, top_anchor[::-1], bottom_anchor[::-1], mode="constant"
        )
        points = _line_profile_coordinates(
            top_anchor[::-1], bottom_anchor[::-1]
        ).swapaxes(0, 1)[:, ::-1, 0]
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
        diagnostics["profiles"] = hv.Overlay([hv.Curve(tp) for tp in padded_profiles])
    if diagnostics is not None:
        lines_plot = hv.Path(
            [[points[0] + 0.5, points[-1] + 0.5] for points in line_points]
        ).options(color="blue")
        top_line_plot = hv.Points([points[0] + 0.5 for points in line_points]).options(
            color="green"
        )
        bottom_line_plot = hv.Points(
            [points[-1] + 0.5 for points in line_points]
        ).options(color="red")
        anchor_points_plot = hv.Points(anchors + 0.5).options(size=3, color="cyan")
        diagnostics["image_with_lines"] = (
            RevImage(img)
            * lines_plot
            * top_line_plot
            * bottom_line_plot
            * anchor_points_plot
        )
    # 0 / 0
    profiles = np.array(padded_profiles)
    stacked_points = np.array(padded_line_points).swapaxes(0, 1)
    return profiles, stacked_points, anchor_idx


def find_trench_ends(
    img,
    angle,
    rhos,
    margin=15,
    profile_quantile=0.95,
    diagnostics=None,
):
    profiles, stacked_points, anchor_idx = get_trench_line_profiles(
        img, angle, rhos, diagnostics=diagnostics
    )
    stacked_profile = silent_nanquantile(profiles, profile_quantile, axis=0)
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
    top_endpoints_shifted = top_endpoints + 0.5
    bottom_endpoints_shifted = bottom_endpoints + 0.5
    if diagnostics is not None:
        trench_plot = hv.Path(
            [
                [top_endpoint, bottom_endpoint]
                for top_endpoint, bottom_endpoint in zip(
                    top_endpoints_shifted, bottom_endpoints_shifted
                )
            ]
        ).options(color="white")
        top_points_plot = hv.Points(top_endpoints_shifted).options(
            size=3, color="green"
        )
        bottom_points_plot = hv.Points(bottom_endpoints_shifted).options(
            size=3, color="red"
        )
        diagnostics["image_with_trenches"] = (
            RevImage(img) * trench_plot * top_points_plot * bottom_points_plot
        )
    trench_index = np.arange(len(mask))[mask]
    df = pd.DataFrame(
        {
            "top_x": top_endpoints[:, 0],
            "top_y": top_endpoints[:, 1],
            "bottom_x": bottom_endpoints[:, 0],
            "bottom_y": bottom_endpoints[:, 1],
        },
        index=trench_index,
    )
    return df
