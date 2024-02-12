import holoviews as hv
import numpy as np
import pandas as pd

from paulssonlab.image_analysis.geometry import get_image_limits
from paulssonlab.image_analysis.trench_detection.end_finding import find_trench_ends
from paulssonlab.image_analysis.trench_detection.hough import find_periodic_lines
from paulssonlab.image_analysis.trench_detection.peaks import find_periodic_peaks
from paulssonlab.image_analysis.ui import RevImage
from paulssonlab.image_analysis.util import getitem_if_not_none


def _get_trench_bbox(top, bottom, width, trench_idx, x_lim, y_lim):
    offset = np.array([width / 2, 0])
    return np.round(
        np.vstack((top[trench_idx] - offset, bottom[trench_idx] + offset))
    ).astype(np.int16)


def _get_trench_bboxes(trenches, width, x_lim, y_lim, **kwargs):
    top_points = np.vstack((trenches["top_x"], trenches["top_y"])).T
    bottom_points = np.vstack((trenches["bottom_x"], trenches["bottom_y"])).T
    return np.hstack(
        [
            _get_trench_bbox(
                top_points, bottom_points, width, trench_idx, x_lim, y_lim, **kwargs
            )[:, np.newaxis]
            for trench_idx in range(len(top_points))
        ]
    )


def get_trench_bboxes(trenches, width, x_lim, y_lim, **kwargs):
    if len(trenches) <= 1:
        upper_left = lower_right = np.full((len(trenches), 2), np.nan)
    else:
        upper_left, lower_right = _get_trench_bboxes(
            trenches, width, x_lim, y_lim, **kwargs
        )
    bboxes = pd.DataFrame(
        {
            "ul_x": upper_left[:, 0],
            "ul_y": upper_left[:, 1],
            "lr_x": lower_right[:, 0],
            "lr_y": lower_right[:, 1],
        },
        index=trenches.index,
    )
    return bboxes


def plot_trenches(trenches_trenches_df, bboxes=True, lines=False, labels=False):
    plots = []
    if lines:
        top_endpoints = (
            np.vstack(
                (
                    trenches_trenches_df["top_x"].values,
                    trenches_trenches_df["top_y"].values,
                )
            ).T
            + 0.5
        )
        bottom_endpoints = (
            np.vstack(
                (
                    trenches_trenches_df["bottom_x"].values,
                    trenches_trenches_df["bottom_y"].values,
                )
            ).T
            + 0.5
        )
        line_plot = hv.Path(
            [
                np.array([top_endpoint, bottom_endpoint])
                for top_endpoint, bottom_endpoint in zip(
                    top_endpoints, bottom_endpoints
                )
            ]
        ).opts(color="white")
        top_points_plot = hv.Points(top_endpoints).opts(size=3, color="green")
        bottom_points_plot = hv.Points(bottom_endpoints).opts(size=3, color="red")
        plots.extend([line_plot, top_points_plot, bottom_points_plot])
    if bboxes:
        bbox_plot = hv.Rectangles(
            (
                trenches_trenches_df["ul_x"],
                trenches_trenches_df["lr_y"] + 1,
                trenches_trenches_df["lr_x"] + 1,
                trenches_trenches_df["ul_y"],
            )
        ).opts(fill_color=None, line_color="red")
        plots.append(bbox_plot)
    if labels:
        # TODO: labels seem to be broken in holoviews/bokeh
        # 1) bokeh JS error
        # 2) bokeh doesn't allow text size to be set in data coördinates (so it scales with zoom level)
        label_plot = hv.Labels(
            (
                trenches_trenches_df["ul_x"],
                trenches_trenches_df["ul_y"],
                trenches_trenches_df.index.values.astype(str),
            )
        ).opts(text_color="white", text_font_size="10pt", xoffset=2, yoffset=2)
        plots.append(label_plot)
    return hv.Overlay(plots)


def find_trenches(
    img,
    angle=None,
    pitch=None,
    width=None,
    width_to_pitch_ratio=None,
    width_to_line_width_ratio=None,
    max_angle=np.deg2rad(10),
    num_angles=400,
    peak_func=find_periodic_peaks,
    end_finding_func=find_trench_ends,
    diagnostics=None,
    join_info=True,
    return_bboxes=True,
):
    if angle is None:
        theta = np.linspace(-max_angle, max_angle, num_angles)
    else:
        theta = [angle]
    angle, rhos, info, line_info = find_periodic_lines(
        img,
        theta=theta,
        pitch=pitch,
        peak_func=peak_func,
        diagnostics=getitem_if_not_none(diagnostics, "find_periodic_lines"),
    )
    if line_info is not None:
        line_info.columns = [f"line_{c}" for c in line_info.columns]
    trenches_df = end_finding_func(
        img,
        angle,
        rhos,
        diagnostics=getitem_if_not_none(diagnostics, "end_finding"),
    )
    if line_info is not None:
        trenches_df = trenches_df.join(line_info, on="trench_line")
    trenches_df.reset_index(inplace=True)
    trenches_df.rename_axis(index="roi", inplace=True)
    if return_bboxes:
        if (
            sum(
                [
                    width is not None,
                    width_to_pitch_ratio is not None,
                    width_to_line_width_ratio is not None,
                ]
            )
            != 1
        ):
            raise ValueError(
                "must specify exactly one of width, width_to_pitch_ratio, width_to_line_width_ratio"
            )
        if width is None:
            if width_to_pitch_ratio is not None:
                if info is not None and "pitch" in info:
                    width = info["pitch"] * width_to_pitch_ratio
                else:
                    raise ValueError(
                        "peak_func must compute pitch if width_to_pitch_ratio is given"
                    )
            elif width_to_line_width_ratio is not None:
                if "line_widths" not in trenches_df.columns:
                    raise ValueError(
                        "peak_func must compute line_widths if width_to_line_width_ratio is given"
                    )
                width = trenches_df["line_widths"].median() * width_to_line_width_ratio
        info = {**(info or {}), "width": width}
        image_limits = get_image_limits(img.shape)
        trench_bboxes = get_trench_bboxes(trenches_df, width, *image_limits)
        if trench_bboxes is not None:
            trenches_df = pd.concat([trenches_df, trench_bboxes], axis=1)
        if diagnostics is not None:
            diagnostics["bboxes"] = RevImage(img) * plot_trenches(trenches_df)
    if join_info:
        if info is not None:
            info_df = pd.DataFrame.from_records([info])
            info_df = info_df.loc[info_df.index.repeat(len(trenches_df))].reset_index(
                drop=True
            )
            trenches_df = trenches_df.join(info_df)
        return trenches_df
    else:
        return trenches_df, info
