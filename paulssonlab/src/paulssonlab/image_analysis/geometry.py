import numpy as np

ROI_COORDINATE_COLUMNS = set(["top", "bottom", "ul", "lr"])


def get_image_limits(shape):
    x_min = y_min = 0
    y_max, x_max = shape  # note order
    # TODO: what convention should we use, should max be inclusive??
    x_max -= 1
    y_max -= 1
    x_lim = (x_min, x_max)
    y_lim = (y_min, y_max)
    return x_lim, y_lim


def bounding_box(points):
    upper_left_x = min(point[0] for point in points)
    upper_left_y = min(point[1] for point in points)
    lower_right_x = max(point[0] for point in points)
    lower_right_y = max(point[1] for point in points)
    return np.array([[upper_left_x, upper_left_y], [lower_right_x, lower_right_y]])


def iter_roi_crops(img, rois, corner=False):
    index = rois.index.values
    ul_x = rois["ul_x"].values
    ul_y = rois["ul_y"].values
    lr_x = rois["lr_x"].values
    lr_y = rois["lr_y"].values
    for i in range(len(index)):
        roi_idx = index[i]
        crop = img[ul_y[i] : lr_y[i] + 1, ul_x[i] : lr_x[i] + 1]
        if corner:
            yield roi_idx, crop, np.array([ul_x[i], ul_y[i]])
        else:
            yield roi_idx, crop


def get_roi_crop(img, rois, roi_idx):
    i = rois.index.get_loc(roi_idx)
    ul_x = rois["ul_x"].values
    ul_y = rois["ul_y"].values
    lr_x = rois["lr_x"].values
    lr_y = rois["lr_y"].values
    return img[ul_y[i] : lr_y[i] + 1, ul_x[i] : lr_x[i] + 1]


def _coordinate_columns(columns):
    cols_x = set([f"{col}_x" for col in ROI_COORDINATE_COLUMNS]) & set(columns)
    cols_y = set([f"{col}_y" for col in ROI_COORDINATE_COLUMNS]) & set(columns)
    return cols_x, cols_y


def filter_rois(rois, image_limits):
    x_lim = image_limits[0]
    y_lim = image_limits[1]
    cols_x, cols_y = _coordinate_columns(rois.columns)
    return rois[
        np.logical_and.reduce([rois[col].between(*x_lim) for col in cols_x])
        & np.logical_and.reduce([rois[col].between(*y_lim) for col in cols_y])
    ]


def shift_rois(rois, shift):
    cols_x, cols_y = _coordinate_columns(rois.columns)
    coords_x = {col: rois[col].values + shift[0] for col in cols_x}
    coords_y = {col: rois[col].values + shift[1] for col in cols_y}
    return rois.assign(**coords_x, **coords_y)
