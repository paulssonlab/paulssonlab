import numpy as np
import pandas as pd
import skimage.morphology
from .cluster import label_for_trenches
from .hough import find_trench_lines
from .refinement import find_trench_ends
from util import getitem_if_not_none
import common

# FROM: https://stackoverflow.com/questions/23815327/numpy-one-liner-for-combining-unequal-length-np-array-to-a-matrixor-2d-array
def stack_jagged(arys, fill=np.nan):
    return np.array(list(zip_longest(*arys, fillvalue=fill))).T


def stack_jagged_points(arys):
    length = max(len(points) for points in arys)
    return np.array(
        [np.pad(points, [(0, length - len(points)), (0, 0)], "edge") for points in arys]
    ).swapaxes(0, 1)


def find_trenches(img, reindex=True, diagnostics=None):
    img_normalized, img_labels, label_index = label_for_trenches(
        img, diagnostics=getitem_if_not_none(diagnostics, "labeling")
    )
    trench_sets = {}
    for label in label_index:
        label_diagnostics = getitem_if_not_none(diagnostics, "label_{}".format(label))
        img_masked = np.where(
            skimage.morphology.binary_dilation(img_labels == label),
            img_normalized,
            np.percentile(img_normalized, 5),
        )
        angle, anchor_rho, rho_min, rho_max, anchor_info = find_trench_lines(
            img_masked,
            diagnostics=getitem_if_not_none(label_diagnostics, "find_trench_lines"),
        )
        trench_sets[label] = find_trench_ends(
            img_masked,
            angle,
            anchor_rho,
            rho_min,
            rho_max,
            diagnostics=getitem_if_not_none(label_diagnostics, "find_trench_ends"),
        )
        if anchor_info is not None:
            anchor_info.columns = [("info", col) for col in anchor_info.columns]
            trench_sets[label] = trench_sets[label].join(anchor_info, how="left")
            if reindex:
                trench_sets[label].reset_index(drop=True, inplace=True)
    trenches_df = pd.concat(trench_sets)
    trenches_df.index.names = ["trench_set", "trench"]
    return trenches_df
