import numpy as np
import dask
from dask import delayed
import dask.array as da
import skimage
from cytoolz import compose, partial
from numbers import Integral
from matriarch_stub import (
    get_nd2_reader,
    get_nd2_frame,
    get_regionprops,
    map_over_labels,
    repeat_apply,
    gaussian_box_approximation,
    hessian_eigenvalues,
    RevImage,
    zarrify,
)
from util import short_circuit_none, none_to_nans, trim_zeros

# TODO: new
def nd2_to_dask(filename, position, channel, rechunk=True):
    nd2 = get_nd2_reader(filename)
    if not isinstance(channel, Integral):
        # compute channel index once so we don't do it per frame in get_nd2_frame
        channel = nd2.metadata["channels"].index(channel)
    frame0 = get_nd2_frame(filename, position, channel, 0)
    _get_nd2_frame = delayed(get_nd2_frame)
    frames = [
        _get_nd2_frame(filename, position, channel, t) for t in range(nd2.sizes["t"])
    ]
    arrays = [
        da.from_delayed(frame, dtype=frame0.dtype, shape=frame0.shape)
        for frame in frames
    ]
    stack = da.stack(arrays, axis=0)
    if rechunk:
        stack = stack.rechunk({0: "auto"})
    return stack


# TODO: new
def nd2_to_futures(client, filename, position, channel):
    nd2 = get_nd2_reader(filename)
    frames = [
        client.submit(get_nd2_frame, filename, position, channel, t)
        for t in range(nd2.sizes["t"])
    ]
    return frames


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
        img_k1, sigmas=np.arange(0.1, 1.5, 0.2)
    )  # , scale_range=(1,3), scale_step=0.5)
    if diagnostics is not None:
        diagnostics["img_k1_frangi"] = RevImage(img_k1_frangi)
    # # TODO: necessary?
    # img_k1_frangi -= img_k1_frangi.min()
    # img_k1_frangi /= img_k1_frangi.max()
    img_k1_frangi_uint = skimage.img_as_uint(img_k1_frangi)
    if diagnostics is not None:
        diagnostics["img_k1_frangi_uint"] = RevImage(img_k1_frangi_uint)
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
    watershed_labels = skimage.morphology.watershed(
        img_k1, clean_seeds, mask=mask, watershed_line=False, compactness=0.01
    )
    watershed_labels = watershed_labels.astype(dtype)
    if diagnostics is not None:
        diagnostics["watershed_labels"] = RevImage(watershed_labels)
        diagnostics["watershed_labels_permuted"] = RevImage(
            permute_labels(watershed_labels)
        )
    return watershed_labels


def process_file(
    col_to_funcs,
    photobleaching_filename,
    segmentation_filename=None,
    flat_fields=None,
    dark_frame=None,
    time_slice=slice(None),
    position_slice=slice(None),
    array_func=zarrify,
    delayed=True,
):
    if delayed is True:
        delayed = dask.delayed
    elif delayed is False:
        delayed = lambda func, **kwargs: func
    if flat_fields is None:
        flat_fields = {}

    def _correct_frame(dark_frame, flat_field, frame):
        if dark_frame is not None:
            frame = frame - dark_frame
            if flat_field is not None:
                flat_field = flat_field - dark_frame
        if flat_field is not None:
            frame = frame / flat_field
        return frame

    def _get_corrected_nd2_frame(dark_frame, flat_field, *args, **kwargs):
        frame = get_nd2_frame(*args, **kwargs)
        if np.any(frame):
            return _correct_frame(dark_frame, flat_field, frame)
        else:
            return None

    _get_corrected_nd2_frame = delayed(_get_corrected_nd2_frame, pure=True)
    regionprops_func = delayed(get_regionprops, pure=True)
    _transpose = partial(short_circuit_none, np.transpose)
    if array_func is not None:
        _none_array_func = partial(short_circuit_none, array_func)
        _array_func = delayed(_none_array_func, pure=True)
        _transpose = compose(_none_array_func, _transpose)
    _none_transpose = compose(_transpose, none_to_nans)
    _none_transpose = delayed(_none_transpose, pure=True)
    _short_circuit_none = delayed(short_circuit_none, pure=True)
    nd2 = get_nd2_reader(photobleaching_filename)
    num_positions = nd2.sizes.get("v", 1)
    num_timepoints = nd2.sizes.get("t", 1)
    photobleaching_channel = nd2.metadata["channels"][0]
    del nd2
    data = {}
    for position in range(num_positions)[position_slice]:
        # TODO: map over dask arrays efficiently
        photobleaching_frames = nd2_to_dask(
            photobleaching_filename, position, photobleaching_channel
        )
        # TODO: correction
        # photobleaching_frames = [_get_corrected_nd2_frame(dark_frame, flat_fields.get(photobleaching_channel), photobleaching_filename, position, 0, t)
        #                             for t in range(num_timepoints)]
        if not segmentation_filename:
            segmentation_filename = photobleaching_filename
        segmentation_channel = get_nd2_reader(segmentation_filename).metadata[
            "channels"
        ][0]
        # segmentation_frame = photobleaching_frames[0]
        segmentation_frame = _get_corrected_nd2_frame(
            dark_frame,
            flat_fields.get(segmentation_channel),
            segmentation_filename,
            position,
            0,
            0,
        )
        segmentation_nd2 = get_nd2_reader(segmentation_filename)
        if segmentation_nd2.metadata["channels"][0] == "BF":
            segmentation_func = compose(segment, invert)
        else:
            segmentation_func = segment
        # TODO: if we zarrify labels, we need to turn back into ndarray before map_over_labels
        # if array_func is not None:
        #     segmentation_func = compose(array_func, segmentation_func)
        segmentation_func = delayed(segmentation_func, pure=True)
        labels = segmentation_func(segmentation_frame)
        # regionprops = regionprops_func(labels, segmentation_frame)
        regionprops = None
        photobleaching_frames = photobleaching_frames[time_slice]
        # traces = {col: _none_transpose([_short_circuit_none(map_over_labels, func, labels, frame)
        #                                                 for frame in photobleaching_frames])
        #               for col, func in col_to_funcs.items()}
        traces = {
            col: _none_transpose(
                [
                    _short_circuit_none(map_over_labels, func, labels, frame)
                    for frame in photobleaching_frames
                ]
            )
            for col, func in col_to_funcs.items()
        }
        if array_func is not None:
            segmentation_frame = _array_func(segmentation_frame)
            labels = _array_func(labels)
        data[position] = {
            "regionprops": regionprops,
            "traces": traces,
            "labels": labels,
            "segmentation_frame": segmentation_frame,
        }
    return data
