import holoviews as hv
import numpy as np
import pandas as pd

from paulssonlab.image_analysis.geometry import get_image_limits, get_trench_bbox
from paulssonlab.image_analysis.trench_detection import find_trenches
from paulssonlab.image_analysis.ui import RevImage
from paulssonlab.image_analysis.util import getitem_if_not_none


def _get_trench_bboxes(trenches, x_lim, y_lim, width_factor=2, **kwargs):
    width = trenches["widths"].median() * width_factor
    top_points = np.vstack((trenches["top_x"], trenches["top_y"])).T
    bottom_points = np.vstack((trenches["bottom_x"], trenches["bottom_y"])).T
    return np.hstack(
        [
            get_trench_bbox(
                top_points, bottom_points, width, trench_idx, x_lim, y_lim, **kwargs
            )[:, np.newaxis]
            for trench_idx in range(len(top_points))
        ]
    )


def get_trench_bboxes(trenches, x_lim, y_lim, **kwargs):
    if len(trenches) <= 1:
        upper_left = lower_right = np.full((len(trenches), 2), np.nan)
    else:
        upper_left, lower_right = _get_trench_bboxes(trenches, x_lim, y_lim, **kwargs)
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


# TODO: integrate this into find_trenches?
def find_trench_bboxes(img, diagnostics=None, **kwargs):
    trenches = find_trenches(
        img, diagnostics=getitem_if_not_none(diagnostics, "find_trenches"), **kwargs
    )
    image_limits = get_image_limits(img.shape)
    trench_bboxes = get_trench_bboxes(trenches, *image_limits)
    if trench_bboxes is not None:
        trenches = pd.concat([trenches, trench_bboxes], axis=1)
    if diagnostics is not None:
        top_endpoints = np.vstack(
            (trenches["top_x"].values, trenches["top_y"].values)
        ).T
        bottom_endpoints = np.vstack(
            (trenches["bottom_x"].values, trenches["bottom_y"].values)
        ).T
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
        bbox_plot = hv.Rectangles(
            (
                trench_bboxes["ul_x"],
                trench_bboxes["lr_y"],
                trench_bboxes["lr_x"],
                trench_bboxes["ul_y"],
            )
        ).opts(fill_color=None, line_color="red")
        label_plot = hv.Labels(
            (trench_bboxes["ul_x"], trench_bboxes["ul_y"], trenches.index.values)
        ).opts(text_color="white", fontsize=10, xoffset=3, yoffset=3)
        diagnostics["bboxes"] = (
            RevImage(img)
            * trench_plot
            * top_points_plot
            * bottom_points_plot
            * bbox_plot
            * label_plot
        )
    return trenches


def iter_crops(img, trenches):
    index = trenches.index.values
    ul_x = trenches["ul_x"].values
    ul_y = trenches["ul_y"].values
    lr_x = trenches["lr_x"].values
    lr_y = trenches["lr_y"].values
    for i in range(len(index)):
        crop = img[ul_y[i] : lr_y[i] + 1, ul_x[i] : lr_x[i] + 1]
        yield i, crop
