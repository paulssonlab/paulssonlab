import numpy as np
import pandas as pd
import skimage.morphology
from .set_finding import binarize_trench_image, find_trench_sets_by_cutting
from .hough import find_trench_lines
from .peaks import find_periodic_peaks
from .refinement import find_trench_ends
from ..util import getitem_if_not_none
from .. import common

# FROM: https://stackoverflow.com/questions/23815327/numpy-one-liner-for-combining-unequal-length-np-array-to-a-matrixor-2d-array
def stack_jagged(arys, fill=np.nan):
    return np.array(list(zip_longest(*arys, fillvalue=fill))).T


def stack_jagged_points(arys):
    length = max(len(points) for points in arys)
    return np.array(
        [np.pad(points, [(0, length - len(points)), (0, 0)], "edge") for points in arys]
    ).swapaxes(0, 1)


def find_trenches(
    img,
    reindex=True,
    setwise=True,  # TODO: set False by default?
    peak_func=find_periodic_peaks,
    set_finding_func=find_trench_sets_by_cutting,
    diagnostics=None,
):
    labeling_diagnostics = getitem_if_not_none(diagnostics, "labeling")
    img_normalized, img_binarized = binarize_trench_image(
        img,
        diagnostics=getitem_if_not_none(labeling_diagnostics, "binarize_trench_image"),
    )
    angle, anchor_rho, rho_min, rho_max, anchor_info = find_trench_lines(
        img_normalized,
        peak_func=peak_func,
        diagnostics=getitem_if_not_none(labeling_diagnostics, "find_trench_lines"),
    )
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
        if setwise:
            angle, anchor_rho, rho_min, rho_max, anchor_info = find_trench_lines(
                img_masked,
                peak_func=peak_func,
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
        trench_sets[label]["trench_set"] = label
        if anchor_info is not None:
            trench_sets[label] = trench_sets[label].join(anchor_info, how="left")
            if reindex:
                # TODO: what is the purpose of reindexing?
                trench_sets[label].reset_index(drop=True, inplace=True)
    trenches_df = pd.concat(trench_sets.values())
    trenches_df.reset_index(drop=True, inplace=True)
    return trenches_df
