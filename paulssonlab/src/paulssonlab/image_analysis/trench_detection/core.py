import holoviews as hv
import numpy as np
import pandas as pd
import skimage.morphology

from paulssonlab.image_analysis.geometry import get_image_limits
from paulssonlab.image_analysis.trench_detection.hough import find_periodic_lines
from paulssonlab.image_analysis.trench_detection.peaks import find_periodic_peaks
from paulssonlab.image_analysis.trench_detection.refinement import find_trench_ends
from paulssonlab.image_analysis.trench_detection.set_finding import (
    binarize_trench_image,
    find_trench_sets_by_cutting,
)
from paulssonlab.image_analysis.ui import RevImage, overlay_inverted_yaxis
from paulssonlab.image_analysis.util import getitem_if_not_none


def _get_trench_bbox(top, bottom, width, trench_idx, x_lim, y_lim):
    half_width = int(np.ceil(width / 2))
    offset = np.array([half_width, 0])
    return np.vstack((top[trench_idx] - offset, bottom[trench_idx] + offset))


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


def plot_trenches(trenches_df, bboxes=True, lines=False, labels=False):
    plots = []
    if lines:
        top_endpoints = np.vstack(
            (trenches_df["top_x"].values, trenches_df["top_y"].values)
        ).T
        bottom_endpoints = np.vstack(
            (trenches_df["bottom_x"].values, trenches_df["bottom_y"].values)
        ).T
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
                trenches_df["ul_x"],
                trenches_df["lr_y"],
                trenches_df["lr_x"],
                trenches_df["ul_y"],
            )
        ).opts(fill_color=None, line_color="red")
        plots.append(bbox_plot)
    if labels:
        # TODO: labels seem to be broken in holoviews/bokeh
        # 1) bokeh JS error
        # 2) bokeh doesn't allow text size to be set in data coördinates (so it scales with zoom level)
        label_plot = hv.Labels(
            (
                trenches_df["ul_x"],
                trenches_df["ul_y"],
                trenches_df.index.values.astype(str),
            )
        ).opts(text_color="white", text_font_size="10pt", xoffset=3, yoffset=3)
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
    set_finding_func=find_trench_sets_by_cutting,
    diagnostics=None,
    join_info=True,
    return_bboxes=True,
):
    labeling_diagnostics = getitem_if_not_none(diagnostics, "labeling")
    img_normalized, img_binarized = binarize_trench_image(
        img,
        diagnostics=getitem_if_not_none(labeling_diagnostics, "binarize_trench_image"),
    )
    if angle is None:
        theta = np.linspace(-max_angle, max_angle, num_angles)
    else:
        theta = [angle]
    angle, anchor_rho, rho_min, rho_max, info, line_info = find_periodic_lines(
        img_normalized,
        theta=theta,
        pitch=pitch,
        peak_func=peak_func,
        diagnostics=getitem_if_not_none(labeling_diagnostics, "find_periodic_lines"),
    )
    if line_info is not None:
        line_info.columns = [f"line_{c}" for c in line_info.columns]
    img_labels, label_index = set_finding_func(
        img_normalized,
        img_binarized,
        angle,
        anchor_rho,
        rho_min,
        rho_max,
        diagnostics=getitem_if_not_none(labeling_diagnostics, "set_finding"),
    )
    trench_sets = {}
    for label in label_index:
        label_diagnostics = getitem_if_not_none(diagnostics, "label_{}".format(label))
        img_masked = np.where(
            skimage.morphology.binary_dilation(img_labels == label),
            img_normalized,
            np.percentile(img_normalized, 5),
        )
        trench_sets[label] = find_trench_ends(
            img_masked,
            angle,
            anchor_rho,
            rho_min,
            rho_max,
            diagnostics=getitem_if_not_none(label_diagnostics, "find_trench_ends"),
        )
        trench_sets[label]["trench_set"] = label
        if line_info is not None:
            trench_sets[label] = trench_sets[label].join(line_info)
    trenches_df = pd.concat(trench_sets.values())
    trenches_df.reset_index(drop=True, inplace=True)
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
            diagnostics["bboxes"] = overlay_inverted_yaxis(
                RevImage(img) * plot_trenches(trenches_df)
            )
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
