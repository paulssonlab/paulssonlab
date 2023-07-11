import holoviews as hv
import numpy as np
from scipy.ndimage import map_coordinates
from skimage._shared.utils import _fix_ndimage_mode, _validate_interpolation_order
from skspatial.objects import Line, LineSegment

from paulssonlab.image_analysis.geometry import get_image_limits
from paulssonlab.image_analysis.trench_detection.geometry import (
    angled_line,
    intersect_line_with_segment,
)
from paulssonlab.image_analysis.ui import RevImage


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
            # need to give coordinates in (y, x) order
            coordinates = _line_profile_coordinates(top[::-1], bottom[::-1])
            line_points = coordinates.swapaxes(0, 1)[:, ::-1, 0]
            line_profile = profile_line(
                img, coordinates=coordinates, mode="constant", cval=np.nan
            )[..., 0]
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
    min_offset = min(o for o in offsets if o is not None)
    max_padded_length = max(
        [
            len(line_profile) + line_offset
            for line_offset, line_profile in zip(offsets, profiles)
            if line_profile is not None
        ]
    )
    padded_profiles = []
    padded_points = []
    for i, (line_profile, line_points, line_offset) in enumerate(
        zip(profiles, points, offsets)
    ):
        if line_profile is None:
            padded_line_profile = np.full(max_padded_length, np.nan)
            padded_line_points = np.full((max_padded_length, 2), np.nan)
        else:
            left_padding = line_offset - min_offset
            right_padding = max_padded_length - left_padding - len(line_profile)
            padded_line_profile = np.pad(
                line_profile,
                [
                    (left_padding, right_padding),
                    *[(0, 0)] * (line_profile.ndim - 1),
                ],
                "constant",
                constant_values=np.nan,
            )
            padded_line_points = np.pad(
                line_points,
                [
                    (left_padding, right_padding),
                    *[(0, 0)] * (line_points.ndim - 1),
                ],
                "edge",
            )
        padded_profiles.append(padded_line_profile)
        padded_points.append(padded_line_points)
    if diagnostics is not None:
        # TODO: make hv.Path??
        diagnostics["profiles"] = hv.Overlay(
            [hv.Curve(line_profile) for line_profile in padded_profiles]
        )
    padded_profiles = np.array(padded_profiles)
    padded_points = np.array(padded_points).swapaxes(0, 1)
    return padded_profiles, padded_points
