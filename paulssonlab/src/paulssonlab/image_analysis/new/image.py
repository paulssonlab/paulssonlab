import numpy as np
import pandas as pd
from paulssonlab.image_analysis.trench_detection import find_trenches
from paulssonlab.image_analysis.geometry import get_image_limits, get_trench_bbox


def _get_trench_bboxes(trenches, x_lim, y_lim, **kwargs):
    top_points = np.vstack((trenches["top_x"], trenches["top_y"])).T
    bottom_points = np.vstack((trenches["bottom_x"], trenches["bottom_y"])).T
    return np.hstack(
        [
            get_trench_bbox(
                top_points, bottom_points, trench_idx, x_lim, y_lim, **kwargs
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
def find_trench_bboxes(img, overlap=0):
    trenches = find_trenches(img)
    image_limits = get_image_limits(img.shape)
    trench_bboxes = get_trench_bboxes(trenches, *image_limits, overlap=overlap)
    if trench_bboxes is not None:
        trenches = pd.concat([trenches, trench_bboxes], axis=1)
    return trenches


def segment_trenches(img, trenches):
    channel_image[ul_y : lr_y + 1, ul_x : lr_x + 1]
