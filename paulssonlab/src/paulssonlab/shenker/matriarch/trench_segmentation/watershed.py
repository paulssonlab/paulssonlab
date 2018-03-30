import numpy as np
import holoviews as hv
import scipy.stats
import skimage
import skimage.morphology
from .core import hessian_eigenvalues
from ui import RevImage
from util import repeat_apply


def _trench_img(img):
    return RevImage(
        img.T
    )  # .options(height=10, width=100)#.options(height=int(img.shape[1]*0.6), width=int(img.shape[0]*0.4))


def uniformize_trench_intensity(img):
    # 0.2 exponent: don't distinguish between medium and high stddev, but distinguish between zero and medium
    std = img.std(axis=-1) ** 0.2
    return img * (std / std.max() / img.mean(axis=-1))[..., :, np.newaxis]


def normalize_trench_intensity(img):
    # 0.2 exponent: don't distinguish between medium and high stddev, but distinguish between zero and medium
    std = img.std(axis=-1) ** 0.2
    normalization = (img.max(axis=-1) * (std / std.sum())).sum()
    return img / normalization


# def segment_trench(img, threshold=_find_holeless_threshold, diagnostics=None):
def segment_trench(img, threshold=4e-2, diagnostics=None):
    if diagnostics is not None:
        diagnostics["img"] = _trench_img(img)
    img = normalize_trench_intensity(uniformize_trench_intensity(img))
    if diagnostics is not None:
        diagnostics["img_normalized"] = _trench_img(img)
    img_scaled = skimage.transform.pyramid_expand(img, upscale=2, multichannel=False)
    if diagnostics is not None:
        diagnostics["img_scaled"] = _trench_img(img_scaled)
    img_blurred = skimage.filters.gaussian(img, 0.5)
    if diagnostics is not None:
        diagnostics["img_blurred"] = _trench_img(img_scaled)
    img_k2 = hessian_eigenvalues(img_blurred)[0]
    # img_k2 = normalize_trench_intensity(uniformize_trench_intensity(img_k2))
    # img_k2 = repeat_apply(skimage.morphology.dilation, 1)(img_k2)
    if diagnostics is not None:
        diagnostics["img_k2"] = _trench_img(img_k2)
    #     img_k2 = skimage.transform.rescale(img_k2, 3)
    #     if diagnostics is not None:
    #         diagnostics['img_k2_rescaled'] = _trench_img(img_k2)
    #     if callable(threshold):
    #         threshold_value = threshold(img_k2)
    #     else:
    #         threshold_value = threshold
    #    bin_img = img_k2 > threshold_value
    bin_img = (-threshold < img_k2) & (img_k2 < threshold)
    bin_img = skimage.morphology.binary_erosion(bin_img)
    bin_img = skimage.morphology.binary_dilation(bin_img)
    # bin_img = skimage.morphology.binary_erosion(bin_img)
    # bin_img = skimage.morphology.binary_dilation(bin_img)
    bin_img = skimage.segmentation.clear_border(bin_img, in_place=True)
    if diagnostics is not None:
        # diagnostics['threshold_value'] = threshold_value
        diagnostics["bin_img"] = _trench_img(bin_img)
    # labeled_img = skimage.morphology.label(bin_img)#(1-bin_img)
    # if diagnostics is not None:
    #     diagnostics['labeled_img'] = _trench_img(labeled_img)
    # if diagnostics is not None:
    #     diagnostics['all_labels'] = _trench_img(all_labels)
    # distances = scipy.ndimage.distance_transform_edt(all_labels)
    # if diagnostics is not None:
    #     diagnostics['distances'] = _trench_img(distances)
    # eroded_distances = repeat_apply(skimage.morphology.erosion, 10)(distances)
    # if diagnostics is not None:
    #     diagnostics['eroded_distances'] = _trench_img(eroded_distances)
    # background = all_labels != (labeled_img != 0)
    # if diagnostics is not None:
    #     diagnostics['background'] = _trench_img(background)
    # labeled_distances = skimage.morphology.label(distances != 0 + background)
    # if diagnostics is not None:
    #     diagnostics['labeled_distances'] = _trench_img(labeled_distances)
    # labeled_distances_clean = skimage.morphology.remove_small_objects(labeled_distances, 20)
    # if diagnostics is not None:
    #     diagnostics['labeled_distances_clean'] = _trench_img(labeled_distances_clean)
    clean_seeds = skimage.morphology.label(
        skimage.morphology.remove_small_objects(bin_img, 20)
    )
    if diagnostics is not None:
        diagnostics["clean_seeds"] = _trench_img(clean_seeds)
    mask = img > skimage.filters.threshold_otsu(img)
    if diagnostics is not None:
        diagnostics["mask"] = _trench_img(mask)
    watershed_labels = skimage.morphology.watershed(
        img_k2, clean_seeds, mask=mask, watershed_line=True
    )
    if diagnostics is not None:
        diagnostics["watershed_labels"] = _trench_img(watershed_labels)
    watershed_labels_eroded = repeat_apply(skimage.morphology.erosion, 1)(
        watershed_labels
    )
    # watershed_labels_rescaled = skimage.transform.rescale(watershed_labels, 1/6, anti_aliasing=False)
    if diagnostics is not None:
        diagnostics["watershed_labels_eroded"] = _trench_img(watershed_labels_eroded)
    return watershed_labels_eroded  # _rescaled
