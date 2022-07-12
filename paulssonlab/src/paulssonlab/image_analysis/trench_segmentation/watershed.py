import numpy as np
import holoviews as hv
import scipy.stats
import skimage
import skimage.morphology
import skimage.segmentation
from ..image import hessian_eigenvalues, permute_labels
from ..ui import RevImage, _trench_img
from ..util import repeat_apply
from .. import common


def uniformize_trench_intensity(img):
    # 0.2 exponent: don't distinguish between medium and high stddev, but distinguish between zero and medium
    std = img.std(axis=-1) ** 0.2
    return img * (std / std.max() / img.mean(axis=-1))[..., :, np.newaxis]


def normalize_trench_intensity(img):
    # 0.2 exponent: don't distinguish between medium and high stddev, but distinguish between zero and medium
    std = img.std(axis=-1) ** 0.2
    normalization = (img.max(axis=-1) * (std / std.sum())).sum()
    return img / normalization


def segment(img, dtype=np.uint16, diagnostics=None):
    img = skimage.img_as_float(img)
    if diagnostics is not None:
        diagnostics["img"] = RevImage(img)
    ##img_scaled = skimage.transform.pyramid_expand(img, upscale=2,
    ##                                              multichannel=False)
    ##if diagnostics is not None:
    ##    diagnostics['img_scaled'] = RevImage(img_scaled)
    img_blurred = skimage.filters.gaussian(img, 0.5)
    if diagnostics is not None:
        diagnostics["img_blurred"] = RevImage(img_blurred)
    mask = img_blurred > skimage.filters.threshold_otsu(img_blurred)
    mask = skimage.morphology.remove_small_objects(mask, 5)
    if diagnostics is not None:
        diagnostics["mask"] = RevImage(mask)
    mask_labels = skimage.morphology.label(mask)
    # img_normalized = normalize_componentwise(img, mask_labels)
    # if diagnostics is not None:
    #    diagnostics['img_normalized'] = RevImage(img_normalized)
    img_k1 = hessian_eigenvalues(img)[0]
    # TODO: necessary?
    img_k1 -= img_k1.min()
    img_k1 /= img_k1.max()
    if diagnostics is not None:
        diagnostics["img_k1"] = RevImage(img_k1)
    # img_k1_frangi = skimage.filters.frangi(img_k1, sigmas=np.arange(0.1,1.5,0.5))#, scale_range=(1,3), scale_step=0.5)
    img_k1_frangi = skimage.filters.frangi(
        img_k1, sigmas=np.arange(0.2, 3, 0.2)
    )  # , scale_range=(1,3), scale_step=0.5)
    if diagnostics is not None:
        diagnostics["img_k1_frangi"] = RevImage(img_k1_frangi)
    # # TODO: necessary?
    img_k1_frangi -= img_k1_frangi.min()
    img_k1_frangi /= img_k1_frangi.max()
    # img_k1_frangi_uint = skimage.img_as_uint(img_k1_frangi)
    # if diagnostics is not None:
    #     diagnostics['img_k1_frangi_uint'] = RevImage(img_k1_frangi_uint)
    # selem = skimage.morphology.disk(20)
    # img_k1_frangi_thresh = skimage.filters.rank.otsu(img_k1_frangi_uint, selem)
    img_k1_frangi_thresh = skimage.filters.threshold_local(
        img_k1_frangi, block_size=23, mode="nearest"
    )
    if diagnostics is not None:
        diagnostics["img_k1_frangi_thresh"] = RevImage(img_k1_frangi_thresh)
    # img_k1_frangi_thresh_blurred = skimage.filters.gaussian(img_k1_frangi_thresh, 0)
    # if diagnostics is not None:
    #     diagnostics['img_k1_frangi_thresh_blurred'] = RevImage(img_k1_frangi_thresh_blurred)
    # img_thresh = img_k1_frangi > img_k1_frangi_thresh_blurred
    img_thresh = img_k1_frangi > img_k1_frangi_thresh
    # img_thresh = img_k1_frangi > skimage.filters.threshold_otsu(img_k1_frangi)
    if diagnostics is not None:
        diagnostics["img_thresh"] = RevImage(img_thresh)
    img_thresh_masked = img_thresh * mask
    if diagnostics is not None:
        diagnostics["img_thresh_masked"] = RevImage(img_thresh_masked)
    # img_thresh_eroded = repeat_apply(skimage.morphology.erosion, 0)(img_thresh_masked)
    # if diagnostics is not None:
    #     diagnostics['img_thresh_eroded'] = RevImage(img_thresh_eroded)
    # clean_seeds = skimage.morphology.label(skimage.morphology.remove_small_objects(img_thresh_eroded, 5))
    clean_seeds = skimage.morphology.label(
        skimage.morphology.remove_small_objects(img_thresh_masked, 5)
    )
    if diagnostics is not None:
        diagnostics["clean_seeds"] = RevImage(clean_seeds)
    watershed_labels = skimage.segmentation.watershed(
        img_k1, clean_seeds, mask=mask, watershed_line=False, compactness=0.01
    )
    watershed_labels = watershed_labels.astype(dtype)
    if diagnostics is not None:
        diagnostics["watershed_labels"] = RevImage(watershed_labels)
        diagnostics["watershed_labels_permuted"] = RevImage(
            permute_labels(watershed_labels)
        )
    return watershed_labels
