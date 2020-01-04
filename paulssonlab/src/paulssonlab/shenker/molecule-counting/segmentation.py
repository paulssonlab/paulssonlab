import numpy as np
import holoviews as hv
import dask
from dask import delayed
from dask.delayed import Delayed
import dask.array as da
import skimage
import scipy.ndimage as ndi
from cytoolz import compose, partial
from numbers import Integral
from matriarch_stub import (
    get_nd2_reader,
    get_nd2_frame,
    get_regionprops,
    repeat_apply,
    gaussian_box_approximation,
    hessian_eigenvalues,
    RevImage,
    zarrify,
)
from util import conditional

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
        stack = stack.rechunk({0: "auto" if rechunk is True else rechunk})
    return stack


def aggregate(func, labels, ary):
    """Applies a reduction func to groups of pixels in ary aggregated by labels.
    labels should be two-dimensional, ary must be at least two-dimensional.
    """
    keys = labels.ravel()
    sorter = np.argsort(keys, kind="mergesort")
    sorted_ = keys[sorter]
    flag = sorted_[:-1] != sorted_[1:]
    slices = np.concatenate(([0], np.flatnonzero(flag) + 1, [keys.size]))
    unique = sorted_[slices[:-1]]
    values = ary.reshape((*ary.shape[:-2], -1))[:, sorter]
    test_eval = func(np.ones((1,) * values.ndim))

    def agg(x):
        groups = np.split(x, slices[1:-1], axis=1)
        reductions = [func(x) for x in groups]
        stack = np.stack(reductions, axis=0)
        return stack

    if isinstance(ary, dask.array.Array):
        chunks = (unique.shape[0], *test_eval.shape[:-1], values.chunks[0])
        new_axis = tuple(range(len(chunks) - 1))
        groups = values.map_blocks(agg, drop_axis=1, new_axis=new_axis, chunks=chunks)
    else:
        groups = [func(x) for x in np.split(values, slices[1:-1], axis=1)]
    return unique, groups


def aggregate_dask(func, labels, ary, dtype=np.float_):
    """Same as ``aggregate`` but labels can be a dask.Delayed object."""
    # TODO: optimized code path for ufuncs?
    values = ary.reshape((*ary.shape[:-2], -1))
    test_eval = func(np.ones((1,) * values.ndim))

    def preprocess(labels):
        keys = labels.ravel()
        sorter = np.argsort(keys, kind="mergesort")
        sorted_ = keys[sorter]
        flag = sorted_[:-1] != sorted_[1:]
        slices = np.concatenate(([0], np.flatnonzero(flag) + 1, [keys.size]))
        unique = sorted_[slices[:-1]]
        return sorter, slices, unique

    if isinstance(labels, Delayed):
        preprocess = delayed(preprocess, nout=3)
    sorter, slices, unique = preprocess(labels)
    if isinstance(labels, Delayed):
        num_labels = np.nan
    else:
        num_labels = unique.shape[0]

    def agg(x, sorter, slices):
        x = x[:, sorter]
        groups = np.split(x, slices[1:-1], axis=1)
        reductions = [func(x) for x in groups]
        stack = np.stack(reductions, axis=0)
        return stack

    if isinstance(ary, dask.array.Array):
        chunks = (num_labels, *test_eval.shape[:-1], values.chunks[0])
        new_axis = tuple(range(len(chunks) - 1))
        groups = values.map_blocks(
            agg,
            sorter,
            slices,
            drop_axis=1,
            new_axis=new_axis,
            chunks=chunks,
            dtype=dtype,
        )
    else:
        groups = [func(x) for x in np.split(values, slices[1:-1], axis=1)]
    return unique, groups


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
    if img is None:
        return None
    # img = skimage.img_as_float(img)
    if diagnostics is not None:
        diagnostics["img"] = RevImage(img)
    ##img_scaled = skimage.transform.pyramid_expand(img, upscale=2,
    ##                                              multichannel=False)
    ##if diagnostics is not None:
    ##    diagnostics['img_scaled'] = RevImage(img_scaled)
    # img_blurred = skimage.filters.gaussian(img, 0.5)
    img_blurred = gaussian_box_approximation(img, 2)
    if diagnostics is not None:
        diagnostics["img_blurred"] = RevImage(img_blurred)
    hist, bin_edges = np.histogram(img_blurred.flat, bins=1024)
    hist_idx = np.argmax(hist)
    threshold_scale_factor = 1.5
    # threshold = bin_edges[idx] * threshold_scale_factor
    # threshold = skimage.filters.threshold_otsu(img_blurred)
    num_thresholds = 50
    max_size = img.size / 10
    thresholds = 10 ** np.linspace(
        np.log10(bin_edges[hist_idx]), np.log10(img_blurred.max()), num_thresholds
    )
    threshold_metrics = []
    for thresh in thresholds:
        labels, _ = ndi.label(img_blurred > thresh)
        sizes = np.bincount(labels.flat)
        # metric = np.median(sizes[min_size:])
        metric = np.median(sizes[sizes < max_size])
        threshold_metrics.append(metric)
    threshold_metrics = np.array(threshold_metrics)
    threshold_metrics[np.isinf(threshold_metrics)] = 0
    threshold_metrics[np.isnan(threshold_metrics)] = 0
    threshold_idx = np.argmax(threshold_metrics)
    threshold = thresholds[threshold_idx]
    if diagnostics is not None:
        # TODO: hv.Histogram doesn't work with logy=True (SEE: https://github.com/holoviz/holoviews/issues/2591)
        # diagnostics['histogram'] = hv.Histogram((bin_edges, hist)) * hv.VLine(threshold).options(color='red',logx=True,logy=True)
        diagnostics["histogram"] = (
            hv.Histogram((bin_edges, np.log10(hist + 1)))
            * hv.VLine(threshold).options(color="red")
            * hv.VLine(bin_edges[hist_idx]).options(color="green")
        ).options(logx=True)
        diagnostics["threshold_scale_factor"] = threshold_scale_factor
        diagnostics["threshold"] = threshold
        diagnostics["threshold_metrics"] = hv.Curve(
            (thresholds, threshold_metrics)
        ) * hv.VLine(threshold).options(color="red")
    mask = img_blurred > threshold
    mask = skimage.morphology.remove_small_objects(mask, 5)
    del img_blurred
    # TODO: is this helpful??
    mask = skimage.morphology.binary_erosion(skimage.morphology.binary_dilation(mask))
    if diagnostics is not None:
        diagnostics["mask"] = RevImage(mask)
    # ndimage.label can work in place, output non-int64 dtypes,
    # and raise an error if the number of labels exceeds the dtype's max value
    # whereas skimage.morphology.label cannot
    mask_labels = np.empty_like(mask, dtype=dtype)
    ndi.label(mask, output=mask_labels)
    # img_normalized = normalize_componentwise(img, mask_labels)
    # if diagnostics is not None:
    #    diagnostics['img_normalized'] = RevImage(img_normalized)
    # TODO: can we compute hessian and frangi without converting to float?
    # for now, we convert to float32 so hessian_eigenvalues doesn't convert it to float64
    img_k1 = hessian_eigenvalues(skimage.img_as_float32(img))[0]
    # TODO: necessary?
    # img_k1 -= img_k1.min()
    # img_k1 /= img_k1.max()
    if diagnostics is not None:
        diagnostics["img_k1"] = RevImage(img_k1)
    # img_k1_frangi = skimage.filters.frangi(img_k1, sigmas=np.arange(0.1,1.5,0.5))#, scale_range=(1,3), scale_step=0.5)
    # img_k1_frangi = skimage.filters.frangi(img_k1, sigmas=np.arange(0.1,1.5,0.2))#, scale_range=(1,3), scale_step=0.5)
    # TODO: frangi breaks if input is float32, why??
    img_k1_frangi = skimage.filters.frangi(
        skimage.img_as_float64(img_k1), sigmas=np.arange(1, 6, 2)
    ).astype(np.float32)
    if diagnostics is not None:
        diagnostics["img_k1_frangi"] = RevImage(img_k1_frangi)
    # # TODO: necessary?
    # img_k1_frangi -= img_k1_frangi.min()
    # img_k1_frangi /= img_k1_frangi.max()
    # img_k1_frangi_uint = skimage.img_as_uint(img_k1_frangi)
    # if diagnostics is not None:
    #     diagnostics['img_k1_frangi_uint'] = RevImage(img_k1_frangi_uint)
    # selem = skimage.morphology.disk(20)
    # img_k1_frangi_thresh = skimage.filters.rank.otsu(img_k1_frangi_uint, selem)
    # img_k1_frangi_thresh = skimage.filters.threshold_local(img_k1_frangi, block_size=23, mode='nearest')
    img_k1_frangi_thresh = np.empty_like(img_k1_frangi)
    ndi.gaussian_filter(img_k1_frangi, 4, output=img_k1_frangi_thresh, mode="nearest")
    if diagnostics is not None:
        diagnostics["img_k1_frangi_thresh"] = RevImage(img_k1_frangi_thresh)
    # img_k1_frangi_thresh_blurred = skimage.filters.gaussian(img_k1_frangi_thresh, 0)
    # if diagnostics is not None:
    #     diagnostics['img_k1_frangi_thresh_blurred'] = RevImage(img_k1_frangi_thresh_blurred)
    # img_thresh = img_k1_frangi > img_k1_frangi_thresh_blurred
    img_thresh = img_k1_frangi > img_k1_frangi_thresh
    del img_k1_frangi, img_k1_frangi_thresh
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
    clean_seeds = np.empty_like(mask, dtype=dtype)
    num_labels = ndi.label(
        skimage.morphology.remove_small_objects(img_thresh_masked, 30),
        output=clean_seeds,
    )
    if diagnostics is not None:
        diagnostics["num_labels"] = num_labels
        diagnostics["clean_seeds"] = RevImage(clean_seeds)
    watershed_labels = skimage.morphology.watershed(
        img_k1, clean_seeds, mask=mask, watershed_line=False, compactness=0.01
    )
    if diagnostics is not None:
        diagnostics["watershed_labels"] = RevImage(watershed_labels)
        diagnostics["watershed_labels_permuted"] = RevImage(
            permute_labels(watershed_labels)
        )
    return watershed_labels


def measure_photobleaching(
    photobleaching_filename,
    position,
    photobleaching_channel,
    labels,
    time_slice=slice(None),
    rechunk=True,
):
    if labels is None:
        return None
    photobleaching_frames = nd2_to_dask(
        photobleaching_filename, position, photobleaching_channel, rechunk=rechunk
    )
    photobleaching_frames = photobleaching_frames[time_slice]
    mean_traces = aggregate(partial(np.mean, axis=-1), labels, photobleaching_frames)[
        1
    ].compute(scheduler="single-threaded")
    traces = {"mean": mean_traces}
    return traces


def cluster_nd2_by_positions(filenames, tol=10, ignored_channels=[]):
    positions = {}
    for filename in filenames:
        nd2 = get_nd2_reader(filename)
        xs = nd2._parser._raw_metadata.x_data
        ys = nd2._parser._raw_metadata.y_data
        for d in (xs, ys):
            if not np.allclose(np.array(d) - d[0], 0):
                raise ValueError(
                    "expected constant x/y stage position in {}".format(filename)
                )
        x = xs[0]
        y = ys[0]
        channels = nd2.metadata["channels"]
        if len(channels) != 1:
            raise ValueError("expected exactly one channel: {}".format(channels))
        if channels[0] in ignored_channels:
            continue
        matched = False
        for pos in positions:
            if np.sqrt((x - pos[0]) ** 2 + (y - pos[1]) ** 2) <= tol:
                if channels[0] in positions[pos]:
                    raise ValueError(
                        "duplicate channel: {} in {}, conflicts with {}".format(
                            channels[0], filename, positions[pos][channels[0]]
                        )
                    )
                else:
                    positions[pos][channels[0]] = filename
                    matched = True
                    break
        if not matched:
            positions[(x, y)] = {channels[0]: filename}
    return positions


def process_photobleaching_file(
    col_to_funcs,
    photobleaching_filename,
    photobleaching_channel=None,
    segmentation_filename=None,
    segmentation_channel=None,
    time_slice=slice(None),
    position_slice=slice(None),
    delayed=True,
    rechunk=True,
    segmentation_frame_filter=None,
    segmentation_labels_filter=None,
):
    if delayed is True:
        delayed = dask.delayed(pure=True)
    elif delayed is False:
        delayed = lambda func, **kwargs: func
    nd2 = get_nd2_reader(photobleaching_filename)
    num_positions = nd2.sizes.get("v", 1)
    if photobleaching_channel is None:
        photobleaching_channels = nd2.metadata["channels"]
        if len(photobleaching_channels) != 1:
            raise ValueError(
                "expected only one photobleaching channel: {} in {}".format(
                    photobleaching_channels, photobleaching_filename
                )
            )
        photobleaching_channel = photobleaching_channels[0]
    if not segmentation_filename:
        segmentation_filename = photobleaching_filename
    if segmentation_channel is None:
        segmentation_channels = get_nd2_reader(segmentation_filename).metadata[
            "channels"
        ]
        if len(segmentation_channels) != 1:
            raise ValueError(
                "expected only one segmentation channel: {} in {}".format(
                    segmentation_channels, segmentation_filename
                )
            )
        segmentation_channel = segmentation_channels[0]
    if segmentation_channel == "BF":
        segmentation_func = compose(segment, invert)
    else:
        segmentation_func = segment
    data = {}
    for position in range(num_positions)[position_slice]:
        segmentation_frame = delayed(get_nd2_frame)(
            segmentation_filename, position, segmentation_channel, 0
        )
        good_segmentation_frame = delayed(segmentation_frame_filter)(segmentation_frame)
        segmentation_frame_filtered = delayed(conditional)(
            good_segmentation_frame, segmentation_frame, None
        )
        labels = delayed(segmentation_func)(segmentation_frame_filtered)
        good_segmentation_labels = delayed(segmentation_labels_filter)(
            labels, segmentation_frame
        )
        labels_filtered = delayed(conditional)(good_segmentation_labels, labels, None)
        traces = delayed(measure_photobleaching)(
            photobleaching_filename,
            position,
            photobleaching_channel,
            labels_filtered,
            time_slice=time_slice,
            rechunk=rechunk,
        )
        regionprops = delayed(get_regionprops)(
            labels_filtered, segmentation_frame_filtered
        )
        data[position] = {
            "traces": traces,
            "segmentation_frame": segmentation_frame,
            "labels": labels,
            "regionprops": regionprops,
        }
    return data
