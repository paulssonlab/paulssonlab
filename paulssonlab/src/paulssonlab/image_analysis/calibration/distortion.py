import holoviews as hv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import skimage
from sklearn.neighbors import KDTree

from paulssonlab.image_analysis.image import scale_image


def segment_puncta(
    img,
    low_sigma=1,
    high_sigma=50,
    expand=2,
):
    img_dog = skimage.filters.difference_of_gaussians(img, low_sigma, high_sigma)
    img_mask = img_dog > skimage.filters.threshold_otsu(img_dog)
    labels = skimage.measure.label(img_mask)
    if expand is not None:
        labels = skimage.segmentation.expand_labels(labels, expand)
    return labels


def measure_puncta(img, labels=None, segment_kwargs={}):
    if labels is None:
        labels = segment_puncta(img, **segment_kwargs)
    background = np.median(img[labels == 0])
    # subtract background (note that some pixels may be negative)
    img_bgsub = img - background
    # normalize image
    img_bgsub /= img_bgsub.max()
    df = pd.DataFrame(
        skimage.measure.regionprops_table(
            labels,
            img_bgsub,
            properties=(
                "label",
                "centroid_weighted",
                "area",
                "moments_weighted_central",
            ),
        )
    )
    # TODO: check for off-by-one
    center_y = (img.shape[0] + 1) / 2
    center_x = (img.shape[1] + 1) / 2
    df["x"] = df["centroid_weighted-1"]
    df["y"] = df["centroid_weighted-0"]
    df["radius"] = np.sqrt((df["x"] - center_x) ** 2 + (df["y"] - center_y) ** 2)
    return df


def df_to_coords(df):
    return np.stack((df["x"], df["y"]), axis=1)


def nearest_neighbors_df(df, df2=None):
    query_coords = df_to_coords(df)
    if df2 is None:
        ref_coords = None
    else:
        ref_coords = df_to_coords(df2)
    return nearest_neighbors(query_coords, ref_coords)


def nearest_neighbors(query_coords, ref_coords=None):
    if ref_coords is None:
        ref_coords = query_coords
        # get nearest non-identical neighbor
        neighbor_idx = 1
    else:
        neighbor_idx = 0
    kdtree = KDTree(ref_coords)
    dists, idxs = kdtree.query(query_coords, k=neighbor_idx + 1)
    dists = dists[:, neighbor_idx]
    idxs = idxs[:, neighbor_idx]
    return dists, idxs


def _filter_puncta_df(
    df,
    min_dist=5,  # 15 is maybe better for PSF-fitting
    max_intensity_factor=1.3,
):
    # filter out points with nearest neighbors closer than min_dist
    dists, idxs = nearest_neighbors_df(df)
    max_intensity = df["moments_weighted_central-0-0"].median() * max_intensity_factor
    bad_mask = (dists < min_dist) | (df["moments_weighted_central-0-0"] > max_intensity)
    bad_labels = df["label"][bad_mask]
    df_filtered = df[~bad_mask]
    return df_filtered, bad_labels


def find_puncta(
    img,
    labels=None,
    filter=True,
    segment_kwargs={},
    filter_kwargs={},
    return_labels=False,
):
    if labels is None:
        labels = segment_puncta(img, **segment_kwargs)
    df = measure_puncta(img, labels)
    if filter:
        df, bad_labels = _filter_puncta_df(df, **filter_kwargs)
        if return_labels:
            labels = labels.copy()
            img_mask = np.isin(labels, bad_labels)
            labels[img_mask] = -labels[img_mask]  # flip sign
    if return_labels:
        return df, labels
    else:
        return df


def translate_df(df, offset):
    return df.assign(x=df["x"] + offset[0], y=df["y"] + offset[1])


def plot_puncta(
    img=None,
    labels=None,
    df=None,
    df2=None,
    scale=True,
    filter_labels=True,
):
    if img is not None:
        img = scale_image(img, scale=scale)
    if labels is None:
        plot_img = img
    else:
        if filter_labels:
            labels = labels.copy()
            labels[labels < 0] = 0
        plot_img = skimage.color.label2rgb(labels, img)
    if img is not None or labels is not None:
        plt.imshow(plot_img, cmap="gray")
    if df is not None:
        if plot_img is not None:
            df = df[
                df["x"].between(0, plot_img.shape[1])
                & df["y"].between(0, plot_img.shape[0])
            ]
        plt.plot(
            df["x"],
            df["y"],
            marker="o",
            mfc="none",
            c="w",
            markersize=4,
            lw=0,
            markeredgewidth=0.5,
            alpha=0.5,
        )
    if df2 is not None:
        if plot_img is not None:
            df2 = df2[
                df2["x"].between(0, plot_img.shape[1])
                & df2["y"].between(0, plot_img.shape[0])
            ]
        plt.plot(
            df2["x"],
            df2["y"],
            marker="x",
            mfc="none",
            c="yellow",
            markersize=2,
            lw=0,
            markeredgewidth=0.5,
            alpha=0.5,
        )


def vectorfield_difference(a, b):
    delta = b - a
    return hv.VectorField.from_uv((a[:, 0], a[:, 1], delta[:, 0], delta[:, 1]))
