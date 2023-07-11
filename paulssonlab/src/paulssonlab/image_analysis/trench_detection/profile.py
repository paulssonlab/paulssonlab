import holoviews as hv
import numpy as np
import pandas as pd
from scipy.ndimage import map_coordinates
from skimage._shared.utils import _fix_ndimage_mode, _validate_interpolation_order
from skspatial.objects import Line, LineSegment

from paulssonlab.image_analysis.geometry import get_image_limits
from paulssonlab.image_analysis.image import hough_bounds
from paulssonlab.image_analysis.misc.holoborodko_diff import holo_diff
from paulssonlab.image_analysis.trench_detection.geometry import (
    angled_line,
    edge_point,
    intersect_line_with_segment,
    trench_anchors,
)
from paulssonlab.image_analysis.ui import RevImage
from paulssonlab.image_analysis.workflow import points_dataframe
from paulssonlab.util.numeric import silent_nanquantile


def _unique_close(objs):
    if not len(objs):
        return objs
    unique = [objs[0]]
    for obj in objs[1:]:
        close = False
        for existing in unique:
            if existing.is_close(obj):
                close = True
                break
        if not close:
            unique.append(obj)
    return unique


# ADAPTED FROM: skimage.measure.profile
def _line_profile_coordinates(src, dst, linewidth=1):
    src_row, src_col = src = np.asarray(src, dtype=float)
    dst_row, dst_col = dst = np.asarray(dst, dtype=float)
    d_row, d_col = dst - src
    theta = np.arctan2(d_row, d_col)

    # we add one above because we include the last point in the profile
    # (in contrast to standard numpy indexing)
    length = int(np.ceil(np.hypot(d_row, d_col) + 1))

    line_row = np.linspace(src_row, dst_row, length)
    line_col = np.linspace(src_col, dst_col, length)

    # we subtract 1 from linewidth to change from pixel-counting
    # (make this line 3 pixels wide) to point distances (the
    # distance between pixel centers)
    row_width = (linewidth - 1) * np.cos(theta) / 2
    col_width = (linewidth - 1) * np.sin(-theta) / 2
    row_offsets = np.linspace(-row_width, row_width, linewidth)
    col_offsets = np.linspace(-col_width, col_width, linewidth)
    perp_rows = line_row[:, np.newaxis] + row_offsets[np.newaxis, :]
    perp_cols = line_col[:, np.newaxis] + col_offsets[np.newaxis, :]
    return np.stack([perp_rows, perp_cols])


# ADAPTED FROM: skimage.measure.profile
def profile_line(
    image,
    src=None,
    dst=None,
    coordinates=None,
    linewidth=1,
    order=None,
    mode="constant",
    cval=np.nan,
):
    order = _validate_interpolation_order(image.dtype, order)
    mode = _fix_ndimage_mode(mode)

    if sum([src is not None and dst is not None, coordinates is not None]) != 1:
        raise ValueError(
            "must specify one of coordinates or both src and dst but not both"
        )
    if coordinates is None:
        coordinates = _line_profile_coordinates(src, dst, linewidth=linewidth)

    if image.ndim == 3:
        pixels = [
            map_coordinates(
                image[..., i],
                coordinates,
                prefilter=order > 1,
                order=order,
                mode=mode,
                cval=cval,
            )
            for i in range(image.shape[2])
        ]
        pixels = np.transpose(np.asarray(pixels), (1, 2, 0))
    else:
        pixels = map_coordinates(
            image, coordinates, prefilter=order > 1, order=order, mode=mode, cval=cval
        )
    # The outputted array with reduce_func=None gives an array where the
    # row values (axis=1) are flipped. Here, we make this consistent.
    pixels = np.flip(pixels, axis=1)

    return pixels


def angled_line_profile_endpoints(angle, rho, x_lim, y_lim):
    # ensure pi/2 <= angle < pi/2
    angle = (angle + np.pi / 2) % np.pi - np.pi / 2
    # the top/bottom labels correspond to the coördinate system
    # where the y-axis is inverted (e.g., what RevImage assumes)
    top_line = LineSegment((x_lim[0], y_lim[0]), (x_lim[1], y_lim[0]))
    bottom_line = LineSegment((x_lim[0], y_lim[1]), (x_lim[1], y_lim[1]))
    left_line = LineSegment((x_lim[0], y_lim[0]), (x_lim[0], y_lim[1]))
    right_line = LineSegment((x_lim[1], y_lim[0]), (x_lim[1], y_lim[1]))
    corner = (x_lim[0], y_lim[0])
    baseline = angled_line(corner, angle + np.pi / 2)
    anchor = baseline.to_point(rho)
    profile_line = angled_line(anchor, angle)
    points = []
    for bounding_line in (top_line, bottom_line, left_line, right_line):
        intersection = intersect_line_with_segment(profile_line, bounding_line)
        if intersection is not None:
            points.append(intersection)
    points = list(
        sorted(
            _unique_close(points),
            key=lambda x: x.distance_point(anchor),
        )
    )
    if len(points) != 2:
        return None, None, None
    top, bottom = points
    if angle <= 0:
        offset_corner = (x_lim[0], y_lim[0])
    else:
        offset_corner = (x_lim[1], y_lim[0])
    offset_baseline = angled_line(offset_corner, angle + np.pi / 2)
    offset = offset_baseline.distance_point(top)
    return top, bottom, offset


def angled_profiles(img, angle, rhos, diagnostics=None):
    x_lim, y_lim = get_image_limits(img.shape)
    profiles = []
    points = []
    offsets = []
    for rho in rhos:
        top, bottom, offset = angled_line_profile_endpoints(angle, rho, x_lim, y_lim)
        if top is None:
            line_profile = None
            line_points = None
            line_offset = None
        else:
            line_points = _line_profile_coordinates(top[::-1], bottom[::-1]).swapaxes(
                0, 1
            )[:, ::-1, 0]
            # TODO: use ndi.map_coordinates
            # need to give coordinates in (y, x) order
            line_profile = profile_line(
                img, top[::-1], bottom[::-1], mode="constant", cval=np.nan
            )
            line_offset = int(np.floor(offset))
        profiles.append(line_profile)
        points.append(line_points)
        offsets.append(line_offset)
    if diagnostics is not None:
        lines_plot = hv.Path(
            [
                [line_points[0] + 0.5, line_points[-1] + 0.5]
                for line_points in points
                if line_points is not None
            ]
        ).options(color="blue")
        top_line_plot = hv.Points(
            [line_points[0] + 0.5 for line_points in points if line_points is not None]
        ).options(color="green")
        bottom_line_plot = hv.Points(
            [line_points[-1] + 0.5 for line_points in points if line_points is not None]
        ).options(color="red")
        diagnostics["image_with_lines"] = (
            RevImage(img) * lines_plot * top_line_plot * bottom_line_plot
        )
    max_offset = max(o for o in offsets if o is not None)
    max_stacked_length = max(
        [
            len(line_profile) + line_offset
            for line_offset, line_profile in zip(offsets, profiles)
            if line_profile is not None
        ]
    )
    padded_profiles = []
    padded_points = []
    for i, (profile, points, offset) in enumerate(zip(profiles, line_points, offsets)):
        if profile is None:
            padded_line_profile = np.full(max_stacked_length, np.nan)
            padded_line_points = np.full((max_stacked_length, 2), np.nan)
        else:
            left_padding = max_offset - offset
            right_padding = max_stacked_length - left_padding - len(profile)
            padded_line_profile = np.pad(
                profile,
                *[(0, 0)] * (profile.ndim - 1),
                "constant",
                constant_values=np.nan,
            )
            padded_line_points = np.pad(
                points,
                [(left_padding, right_padding), *[(0, 0)] * (points.ndim - 1)],
                "edge",
            )
        padded_profiles.append(padded_line_profile)
        padded_points.append(padded_line_points)
    if diagnostics is not None:
        # TODO: make hv.Path??
        diagnostics["profiles"] = hv.Overlay(
            [hv.Curve(line_profile) for line_profile in padded_profiles]
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
        coordinates = _line_profile_coordinates(top_anchor[::-1], bottom_anchor[::-1])
        points = coordinates.swapaxes(0, 1)[:, ::-1, 0]
        profile = profile_line(img, coordinates=coordinates, mode="constant")[..., 0]
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
    max_stacked_length = max(
        [len(profile) - offset for offset, profile in zip(offsets, profiles)]
    )
    padded_profiles = []
    padded_line_points = []
    for i, (profile, points, offset) in enumerate(zip(profiles, line_points, offsets)):
        left_padding = offset - min_offset
        right_padding = max_stacked_length - left_padding - len(profile)
        padded_profile = np.pad(
            profile,
            [(left_padding, right_padding), *[(0, 0)] * (profile.ndim - 1)],
            "constant",
            constant_values=np.nan,
        )
        padded_points = np.pad(
            points,
            [(left_padding, right_padding), *[(0, 0)] * (points.ndim - 1)],
            "edge",
        )
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
        diagnostics["image_with_lines"] = (
            RevImage(img) * lines_plot * top_line_plot * bottom_line_plot
        )
    profiles = np.array(padded_profiles)
    stacked_points = np.array(padded_line_points).swapaxes(0, 1)
    return profiles, stacked_points


def find_trench_ends(
    img,
    angle,
    rhos,
    margin=15,
    profile_quantile=0.95,
    diagnostics=None,
):
    profiles, stacked_points = get_trench_line_profiles(
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
