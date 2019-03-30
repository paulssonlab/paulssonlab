import numpy as np
import dask
from dask import delayed
import dask.array as da
import skimage
from cytoolz import compose
from matriarch_stub import (
    get_nd2_reader,
    get_nd2_frame,
    get_regionprops,
    map_over_labels,
    repeat_apply,
    gaussian_box_approximation,
    hessian_eigenvalues,
    RevImage,
)
from util import short_circuit_none, none_to_nans, trim_zeros

# TODO: new
def nd2_to_dask(filename, position, channel):
    nd2 = get_nd2_reader(filename)
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
def segment(img, diagnostics=None):
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
    img_k2 = hessian_eigenvalues(img)[1]
    if diagnostics is not None:
        diagnostics["img_k2"] = RevImage(img_k2)
    # img_k2_frangi = skimage.filters.frangi(img_k2, sigmas=np.arange(0.1,1.5,0.5))#, scale_range=(1,3), scale_step=0.5)
    img_k2_frangi = skimage.filters.frangi(
        img_k2, sigmas=np.arange(0.1, 1.5, 0.2)
    )  # , scale_range=(1,3), scale_step=0.5)
    if diagnostics is not None:
        diagnostics["img_k2_frangi"] = RevImage(img_k2_frangi)
    img_thresh = img_k2_frangi > skimage.filters.threshold_otsu(img_k2_frangi)
    if diagnostics is not None:
        diagnostics["img_thresh"] = RevImage(img_thresh)
    clean_seeds = skimage.morphology.label(
        skimage.morphology.remove_small_objects(img_thresh * mask, 5)
    )
    if diagnostics is not None:
        diagnostics["clean_seeds"] = RevImage(clean_seeds)
    watershed_labels = skimage.morphology.watershed(
        img_k2, clean_seeds, mask=mask, watershed_line=False, compactness=0.01
    )
    watershed_labels = watershed_labels  # .astype(np.uint8)
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
    initial_filename=None,
    final_filename=None,
    flat_fields=None,
    dark_frame=None,
    time_slice=slice(None),
    position_slice=slice(None),
):
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
    _none_transpose = delayed(compose(np.transpose, none_to_nans), pure=True)
    _short_circuit_none = delayed(short_circuit_none, pure=True)
    nd2 = get_nd2_reader(photobleaching_filename)
    num_positions = nd2.sizes.get("v", 1)
    num_timepoints = nd2.sizes.get("t", 1)
    photobleaching_channel = nd2.metadata["channels"][0]
    del nd2
    data = {}
    for position in range(num_positions)[position_slice]:
        # TODO: map over dask arrays efficiently
        # photobleaching_frames = nd2_to_dask(photobleaching_filename, position, 0)
        photobleaching_frames = [
            _get_corrected_nd2_frame(
                dark_frame,
                flat_fields.get(photobleaching_channel),
                photobleaching_filename,
                0,
                0,
                t,
            )
            for t in range(num_timepoints)
        ]
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
        segmentation_func = delayed(segmentation_func, pure=True)
        labels = segmentation_func(segmentation_frame)
        sandwich_frames = []
        if initial_filename is not None:
            sandwich_channel = get_nd2_reader(initial_filename).metadata["channels"][0]
            initial_frame = _get_corrected_nd2_frame(
                dark_frame,
                flat_fields.get(sandwich_channel),
                initial_filename,
                position,
                0,
                0,
            )
            sandwich_frames.append(initial_frame)
            regionprops_frame = initial_frame
        else:
            regionprops_frame = segmentation_frame
        if final_filename is not None:
            sandwich_frames.append(
                _get_corrected_nd2_frame(
                    dark_frame,
                    flat_fields.get(sandwich_channel),
                    final_filename,
                    position,
                    0,
                    0,
                )
            )
        sandwich_traces = {
            col: _none_transpose(
                [
                    _short_circuit_none(map_over_labels, func, labels, frame)
                    for frame in sandwich_frames
                ]
            )
            for col, func in col_to_funcs.items()
        }
        regionprops = regionprops_func(labels, regionprops_frame)
        photobleaching_frames = photobleaching_frames[time_slice]
        traces = {
            col: _none_transpose(
                [
                    _short_circuit_none(map_over_labels, func, labels, frame)
                    for frame in photobleaching_frames
                ]
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
