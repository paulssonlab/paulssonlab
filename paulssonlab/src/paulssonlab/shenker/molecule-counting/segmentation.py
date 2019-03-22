import numpy as np
import dask
import dask.array as da
import skimage
from cytoolz import compose
from matriarch_stub import (
    get_nd2_reader,
    get_nd2_frame,
    get_regionprops,
    nd2_to_dask,
    nd2_to_futures,
    map_over_labels,
    repeat_apply,
    gaussian_box_approximation,
    hessian_eigenvalues,
    RevImage,
)
from util import short_circuit_none, none_to_nans

# FROM: http://emmanuelle.github.io/a-tutorial-on-segmentation.html
def permute_labels(labels):
    label_map = np.concatenate(((0,), np.random.permutation(labels.max()) + 1))
    return label_map[labels]


def invert(img):
    img2 = img.astype(np.float32)
    img2 = np.median(img2) - img2
    img2[img2 < 0] = 0
    return img2


def segment_otsuonly(img):
    img_highpass = img - gaussian_box_approximation(img, 5)
    mask = img > skimage.filters.threshold_otsu(img)
    mask = skimage.morphology.remove_small_objects(mask, 3)
    labels = skimage.measure.label(mask)
    return labels


# TODO: modified
def segment(img, diagnostics=None):
    img = img.astype(
        np.float32
    )  # TODO: automatically cast to float32 if not already float32/64
    if diagnostics is not None:
        diagnostics["img"] = RevImage(img)
    # img = normalize_trench_intensity(uniformize_trench_intensity(img))
    # if diagnostics is not None:
    #    diagnostics['img_normalized'] = RevImage(img)
    ##img_scaled = skimage.transform.pyramid_expand(img, upscale=2,
    ##                                              multichannel=False)
    ##if diagnostics is not None:
    ##    diagnostics['img_scaled'] = RevImage(img_scaled)
    img_blurred = skimage.filters.gaussian(img, 0.5)
    if diagnostics is not None:
        diagnostics["img_blurred"] = RevImage(img_blurred)
    img_k1 = hessian_eigenvalues(img_blurred)[0]
    mask = img_blurred > skimage.filters.threshold_otsu(img_blurred)
    if diagnostics is not None:
        diagnostics["mask"] = RevImage(mask)
    if diagnostics is not None:
        diagnostics["img_k1"] = RevImage(img_k1)
    img_k1_frangi = skimage.filters.frangi(
        img_k1, sigmas=np.arange(0.1, 1.5, 0.5)
    )  # , scale_range=(1,3), scale_step=0.5)
    if diagnostics is not None:
        diagnostics["img_k1_frangi"] = RevImage(img_k1_frangi)
    img_thresh = img_k1_frangi > skimage.filters.threshold_otsu(img_k1_frangi)
    if diagnostics is not None:
        diagnostics["img_thresh"] = RevImage(img_thresh)
    clean_seeds = skimage.morphology.label(
        skimage.morphology.remove_small_objects(img_thresh * mask, 5)
    )
    if diagnostics is not None:
        diagnostics["clean_seeds"] = RevImage(clean_seeds)
    watershed_labels = skimage.morphology.watershed(
        img_k1, clean_seeds, mask=mask, watershed_line=False, compactness=0.01
    )
    watershed_labels = watershed_labels  # .astype(np.uint8)
    if diagnostics is not None:
        diagnostics["watershed_labels"] = RevImage(watershed_labels)
    return watershed_labels


# def segment(img):
#     img_frangi = skimage.filters.frangi(img, scale_range=(0.1,1.5), scale_step=0.1)
#     del img
#     mask = img_frangi < np.percentile(img_frangi, 90)
#     del img_frangi
#     mask = skimage.segmentation.clear_border(mask)
#     mask = repeat_apply(skimage.morphology.erosion, 2)(mask)
#     mask = skimage.morphology.remove_small_objects(mask, 5)
#     labels = skimage.measure.label(mask)
#     # TODO: not even watershedding??
#     return labels


def process_file(
    client,
    col_to_funcs,
    photobleaching_filename,
    segmentation_filename=None,
    initial_filename=None,
    final_filename=None,
    time_slice=slice(None),
):
    nd2 = get_nd2_reader(photobleaching_filename)
    num_positions = nd2.sizes.get("v", 1)
    del nd2
    data = {}
    for position in range(num_positions):
        # TODO: why won't dask array work?
        # photobleaching_frames = nd2_to_dask(photobleaching_filename, position, 0)
        photobleaching_frames = nd2_to_futures(
            client, photobleaching_filename, position, 0
        )
        if segmentation_filename:
            segmentation_frame = client.submit(
                get_nd2_frame, segmentation_filename, position, 0, 0
            )
        else:
            segmentation_filename = photobleaching_filename
            segmentation_frame = photobleaching_frames[0]
        segmentation_nd2 = get_nd2_reader(segmentation_filename)
        if segmentation_nd2.metadata["channels"][0] == "BF":
            segmentation_func = compose(segment, invert)
        else:
            segmentation_func = segment
        # TODO: assuming segmenting in phase
        labels = client.submit(segmentation_func, segmentation_frame)
        sandwich_frames = []
        if initial_filename is not None:
            initial_frame = client.submit(
                get_nd2_frame, initial_filename, position, 0, 0
            )
            sandwich_frames.append(initial_frame)
            regionprops_frame = initial_frame
        else:
            regionprops_frame = segmentation_frame
        if final_filename is not None:
            sandwich_frames.append(
                client.submit(get_nd2_frame, final_filename, position, 0, 0)
            )
        sandwich_traces = {
            col: client.submit(
                compose(np.transpose, none_to_nans),
                [
                    client.submit(
                        short_circuit_none, map_over_labels, func, labels, frame
                    )
                    for frame in sandwich_frames
                ],
            )
            for col, func in col_to_funcs.items()
        }
        regionprops = client.submit(get_regionprops, labels, regionprops_frame)
        # traces = {col: client.submit(np.transpose, client.map(partial(map_over_labels, func, labels),
        #                                                      photobleaching_frames[time_slice]))
        traces = {
            col: client.submit(
                compose(np.transpose, none_to_nans),
                [
                    client.submit(
                        short_circuit_none, map_over_labels, func, labels, frame
                    )
                    for frame in photobleaching_frames[time_slice]
                ],
            )
            for col, func in col_to_funcs.items()
        }
        data[position] = {
            "regionprops": regionprops,
            "sandwich_traces": sandwich_traces,
            "traces": traces,
            "labels": labels,
            "segmentation_frame": segmentation_frame,
        }
    return data
