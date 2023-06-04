import numpy as np
import pandas as pd
import skimage.morphology

from paulssonlab.image_analysis.trench_detection.hough import find_periodic_lines
from paulssonlab.image_analysis.trench_detection.peaks import find_periodic_peaks
from paulssonlab.image_analysis.trench_detection.refinement import find_trench_ends
from paulssonlab.image_analysis.trench_detection.set_finding import (
    binarize_trench_image,
    find_trench_sets_by_cutting,
)
from paulssonlab.image_analysis.util import getitem_if_not_none


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
    max_angle=np.deg2rad(10),
    num_angles=400,
    peak_func=find_periodic_peaks,
    set_finding_func=find_trench_sets_by_cutting,
    diagnostics=None,
    join_info=True,
):
    labeling_diagnostics = getitem_if_not_none(diagnostics, "labeling")
    img_normalized, img_binarized = binarize_trench_image(
        img,
        diagnostics=getitem_if_not_none(labeling_diagnostics, "binarize_trench_image"),
    )
    angle, anchor_rho, rho_min, rho_max, info, line_info = find_periodic_lines(
        img_normalized,
        theta=np.linspace(-max_angle, max_angle, num_angles),
        peak_func=peak_func,
        diagnostics=getitem_if_not_none(labeling_diagnostics, "find_periodic_lines"),
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
