import itertools as it

import holoviews as hv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import skimage
from scipy.optimize import minimize
from sklearn.neighbors import KDTree
from tqdm.auto import tqdm, trange

from paulssonlab.image_analysis.image import radial_distortion, scale_image


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


def optimize_radial_distortion_correction(nd2, **kwargs):
    # def bounds_func(shape):
    #     return ((0, 1e-8), (0, nd2.sizes["x"] - 1), (0, nd2.sizes["y"] - 1))

    def correction_func(coords, params):
        return radial_distortion(
            coords, params[0], (params[1], params[2]), (params[1], params[2])
        )

    # return optimize_correction(nd2, correction_func, initial_params, bounds, **kwargs)
    initial_params = (1e-9, nd2.sizes["x"] / 2, nd2.sizes["y"] / 2)
    bounds = ((0, 1e-8), (0, nd2.sizes["x"] - 1), (0, nd2.sizes["y"] - 1))
    res = optimize_correction(nd2, correction_func, initial_params, bounds, **kwargs)
    return res  # TODO


def _prepare_optimize_correction(
    nd2,
    fov_pairs,
    transform_type=skimage.transform.EuclideanTransform,
    t=0,
    z_idx=0,
    channel_idx=0,
    objective_magnification=20,
    pixel_size=4.25,
):
    xs = np.asarray(nd2._parser._raw_metadata.x_data)
    ys = np.asarray(nd2._parser._raw_metadata.y_data)
    coords_distorted = {}
    transforms = {}
    for fov_pair in fov_pairs:
        fov1, fov2 = fov_pair
        if (coords1 := coords_distorted.get(fov1)) is None:
            img = nd2.get_frame_2D(v=fov1, t=t, z=z_idx, c=channel_idx)
            coords1 = df_to_coords(find_puncta(img))
            del img
            coords_distorted[fov1] = coords1
        if (coords2 := coords_distorted.get(fov2)) is None:
            img = nd2.get_frame_2D(v=fov2, t=t, z=z_idx, c=channel_idx)
            coords2 = df_to_coords(find_puncta(img))
            del img
            coords_distorted[fov2] = coords2
        stage_offset = np.array(
            [
                xs[fov2] - xs[fov1],
                ys[fov2] - ys[fov1],
            ]
        )
        translation_guess = stage_offset * objective_magnification / pixel_size
        # TODO: try TranslationTransform?
        transform = transform_type(translation=translation_guess)
        transforms[fov_pair] = transform
    return coords_distorted, transforms


def optimize_correction(
    nd2,
    correction_func,
    initial_params,
    bounds,
    transform_type=skimage.transform.EuclideanTransform,
    t=0,
    z_idx=0,
    channel_idx=0,
    fov_pairs=None,
    objective_magnification=20,
    pixel_size=4.25,
    max_correspondence_dist=5,
    num_iters=1,
    progress_bar=tqdm,
    prepared=None,
    **kwargs,
):
    if fov_pairs is None:
        fov_pairs = it.pairwise(range(nd2.sizes["v"]))
    fov_pairs = list(fov_pairs)
    if progress_bar is not None:
        fov_pairs_progress = progress_bar(fov_pairs)
    else:
        fov_pairs_progress = fov_pairs
    if prepared is None:
        coords_distorted, transforms = _prepare_optimize_correction(
            nd2,
            fov_pairs_progress,
            transform_type=transform_type,
            t=t,
            z_idx=z_idx,
            channel_idx=channel_idx,
            objective_magnification=objective_magnification,
            pixel_size=pixel_size,
        )
    else:
        coords_distorted, transforms = prepared

    def objective_func(params, correction_func, coords1_all, coords2_all, transform):
        se = 0
        for coords1, coords2 in zip(coords1_all, coords2_all):
            coords1_corrected = correction_func(coords1, params)
            coords2_corrected = correction_func(coords2, params)
            transform.estimate(coords1_corrected, coords2_corrected)
            se += (
                (coords1_corrected - transform.inverse(coords2_corrected)) ** 2
            ).sum()
        rmse = np.sqrt(se / len(coords1_all))
        return rmse

    params = initial_params
    for iter_num in trange(num_iters + 1):
        coords1_correspondences = []
        coords2_correspondences = []
        for fov_pair in fov_pairs:
            fov1, fov2 = fov_pair
            coords1_distorted = coords_distorted[fov1]
            coords2_distorted = coords_distorted[fov2]
            coords1 = correction_func(coords1_distorted, params)
            coords2 = correction_func(coords2_distorted, params)
            transform = transforms[fov_pair]
            coords2_transformed = transform.inverse(coords2)
            correspondence_dists, correspondence_idxs = nearest_neighbors(
                coords1, coords2_transformed
            )
            correspondence_mask = correspondence_dists < max_correspondence_dist
            coords1_correspondence = coords1_distorted[correspondence_mask]
            coords2_correspondence = coords2_distorted[correspondence_idxs][
                correspondence_mask
            ]
            transform.estimate(coords1_correspondence, coords2_correspondence)
            coords1_correspondences.append(coords1_correspondence)
            coords2_correspondences.append(coords2_correspondence)
        if iter_num == num_iters:
            # last iteration: re-estimate transforms/correspondences before returning
            # (this only matters if we return the transforms or correspondences)
            return res
        res = minimize(
            objective_func,
            initial_params,
            **{
                "method": "L-BFGS-B",
                "bounds": bounds,
                "args": (
                    correction_func,
                    coords1_correspondences,
                    coords2_correspondences,
                    transform_type(),
                ),
                **kwargs,
            },
        )
        if not res.success:
            raise Exception(
                f"optimizer not successful (status {res.status}): {res.message}"
            )
        params = res.x


def plot_puncta(
    img=None,
    labels=None,
    coords=None,
    coords2=None,
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
    if coords is not None:
        if isinstance(coords, pd.DataFrame):
            coords_x = coords["x"]
            coords_y = coords["y"]
        else:
            coords_x = coords[:, 0]
            coords_y = coords[:, 1]
        if plot_img is not None:
            coords_mask = (
                (-0.5 <= coords_x)
                & (coords_x <= plot_img.shape[1] - 0.5)
                & (-0.5 <= coords_y)
                & (coords_y <= plot_img.shape[0] - 0.5)
            )
            coords_x = coords_x[coords_mask]
            coords_y = coords_y[coords_mask]
        plt.plot(
            coords_x,
            coords_y,
            marker="o",
            mfc="none",
            c="w",
            markersize=4,
            lw=0,
            markeredgewidth=0.5,
            alpha=0.5,
        )
    if coords2 is not None:
        if isinstance(coords2, pd.DataFrame):
            coords2_x = coords2["x"]
            coords2_y = coords2["y"]
        else:
            coords2_x = coords2[:, 0]
            coords2_y = coords2[:, 1]
        if plot_img is not None:
            coords2_mask = (
                (-0.5 <= coords2_x)
                & (coords2_x <= plot_img.shape[1] - 0.5)
                & (-0.5 <= coords2_y)
                & (coords2_y <= plot_img.shape[0] - 0.5)
            )
            coords2_x = coords2_x[coords2_mask]
            coords2_y = coords2_y[coords2_mask]
        plt.plot(
            coords2_x,
            coords2_y,
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
