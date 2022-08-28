import numpy as np
from matplotlib.colors import hex2color
import dask
import scipy.ndimage as ndi
import skimage
from skimage.transform import SimilarityTransform, warp
import nd2reader
import av
import itertools as it
from numbers import Number
from cytoolz import partial
from tqdm.auto import tqdm
from paulssonlab.image_analysis.workflow import (
    parse_nd2_metadata,
    get_position_metadata,
    get_filename_image_limits,
    get_nd2_frame,
)
from paulssonlab.image_analysis.image import gaussian_box_blur

# TODO: move to paulssonlab.util
def get_delayed(delayed):
    if delayed is True:
        return dask.delayed(pure=True)
    elif delayed is False:
        return lambda func, **kwargs: func
    else:
        return delayed


def composite_channels(imgs, hexcolors, scale=True):
    colors = [hex2color(hexcolor) for hexcolor in hexcolors]
    return _composite_channels(imgs, colors, scale=scale)


def _composite_channels(channel_imgs, colors, scale=True):
    if len(channel_imgs) != len(colors):
        raise ValueError("expecting equal numbers of channels and colors")
    num_channels = len(channel_imgs)
    if scale:
        scaled_imgs = [
            channel_imgs[i] / np.percentile(channel_imgs[i], 99.9)
            for i in range(num_channels)
        ]
        for scaled_img in scaled_imgs:
            np.clip(scaled_img, 0, 1, scaled_img)  # clip in place
    else:
        scaled_imgs = channel_imgs
    imgs_to_combine = [
        scaled_imgs[i][:, :, np.newaxis] * np.array(colors[i])
        for i in range(num_channels)
    ]
    if not len(imgs_to_combine):
        imgs_to_combine = [np.ones(colored_imgs[0].shape)]  # white placeholder
    img = imgs_to_combine[0]
    for img2 in imgs_to_combine[1:]:
        img = 1 - (1 - img) * (1 - img2)
    return img


def colorized_frame(
    get_frame_func,
    filename,
    channels,
    channel_to_color,
    *,
    t=0,
    v=0,
    scaling_funcs=None,
):
    if channels is None:
        raise ValueError("must specify channels")
    imgs = [get_frame_func(filename, v, channel, t) for channel in channels]
    if scaling_funcs:
        for idx in range(len(channels)):
            channel = channels[idx]
            if channel not in scaling_funcs:
                raise ValueError(f"missing scaling_func for {channel}")
            imgs[idx] = scaling_funcs[channel](imgs[idx])
    img = composite_channels(
        imgs,
        [channel_to_color[channel] for channel in channels],
        scale=(not scaling_funcs),
    )
    return img


def rectangles_intersect(ul1, lr1, ul2, lr2):
    return not (
        (ul1[0] > lr2[0]) or (lr1[0] < ul2[0]) or (ul1[1] > lr2[1]) or (lr1[1] < ul2[1])
    )


def scale_around_center(scale, center):
    x, y = center
    return (
        SimilarityTransform(translation=(-x, -y))
        + SimilarityTransform(scale=scale)
        + SimilarityTransform(translation=(x, y))
    )


def output_transformation(input_width, input_height, output_width, output_height):
    width_ratio = input_width / output_width
    height_ratio = input_height / output_height
    scale = max(width_ratio, height_ratio)
    x = -(output_width - input_width / scale) / 2
    y = -(output_height - input_height / scale) / 2
    return SimilarityTransform(translation=(x, y)) + SimilarityTransform(scale=scale)


def mosaic_frame(
    get_frame_func,
    positions,
    image_dims,
    *,
    timepoint=None,
    center=None,
    offset=None,
    scale=1,
    output_dims=(1024, 1024),
    anti_aliasing=True,
    delayed=False,
):
    delayed = get_delayed(delayed)
    if center is not None and offset is not None:
        raise ValueError("can only specify at most one of center, offset")
    if center is None:
        columns = positions["x_idx"].max() - positions["x_idx"].min() + 1
        rows = positions["y_idx"].max() - positions["y_idx"].min() + 1
        center = np.array([image_dims[0] * columns / 2, image_dims[1] * rows / 2])
    if offset is not None:
        center += np.array(offset)
    viewport_transform = output_transformation(*image_dims, *output_dims)
    output_img = delayed(np.zeros)((output_dims[1], output_dims[0], 3))
    viewport_ul = (0, 0)
    viewport_lr = (output_dims[0] - 1, output_dims[1] - 1)  # TODO: off-by-one?
    for (filename, pos_num), position in positions.iterrows():
        frame_corner = (
            -image_dims[0] * position["x_idx"],
            -image_dims[1] * position["y_idx"],
        )
        frame_transform = (
            output_transformation(*image_dims, *output_dims)
            + scale_around_center(1 / scale, (image_dims[0] / 2, image_dims[1] / 2))
            + SimilarityTransform(
                translation=(
                    center[0] - image_dims[0] / 2,
                    center[1] - image_dims[1] / 2,
                )
            )
            + SimilarityTransform(translation=frame_corner)
        )
        frame_ul = frame_transform.inverse((0, 0))[0]
        frame_lr = frame_transform.inverse((image_dims[0] - 1, image_dims[1] - 1))[0]
        visible = rectangles_intersect(viewport_ul, viewport_lr, frame_ul, frame_lr)
        if visible:
            img = delayed(get_frame_func)(t=timepoint, v=pos_num)
            if anti_aliasing:
                # taken from skimage.transform.resize
                anti_aliasing_sigma = (
                    max(0, (1 / scale - 1) / 2) * 2
                )  # TODO: fix fudge factor
                img = delayed(ndi.gaussian_filter)(
                    img, (anti_aliasing_sigma, anti_aliasing_sigma, 0), mode="constant"
                )
            output_img += delayed(warp)(
                img, frame_transform, output_shape=output_dims[::-1], order=3
            )
    return output_img


def mosaic_animate_scale(
    filename,
    scale=1,
    timepoints=None,
    center=None,
    offset=None,
    width=1024,
    height=1024,
    # frame_rate=1, #TODO
    channels=None,
    channel_to_color=None,
    scaling_funcs=None,
    delayed=True,
    progress_bar=tqdm,
    # ignore_exceptions=True,
):
    if channels is None:
        raise ValueError("must specify channels")
    if channel_to_color is None:
        raise ValueError("must specify channel_to_color")
    delayed = get_delayed(delayed)
    # TODO
    # if ignore_exceptions:
    #     excepts_get_nd2_frame = excepts(Exception, get_nd2_frame)
    #     excepts_segmentation_func = excepts(Exception, segmentation_func)
    #     excepts_measure = excepts(Exception, measure)
    # else:
    #     excepts_get_nd2_frame = get_nd2_frame
    #     excepts_segmentation_func = segmentation_func
    #     excepts_measure = measure
    nd2 = nd2reader.ND2Reader(filename)
    nd2s = {filename: nd2 for filename in (filename,)}
    metadata = {
        nd2_filename: parse_nd2_metadata(nd2) for nd2_filename, nd2 in nd2s.items()
    }
    positions = get_position_metadata(metadata)
    # TODO
    # small_positions = positions[(positions["y_idx"] < 3) & (positions["x_idx"] < 3)]
    image_limits = get_filename_image_limits(metadata)
    get_frame_func = partial(
        colorized_frame,
        get_nd2_frame,
        filename,
        channels,
        channel_to_color,
        scaling_funcs=scaling_funcs,
    )
    input_dims = (
        image_limits[filename][0][1] + 1,
        image_limits[filename][1][1] + 1,
    )
    if isinstance(scale, Number):
        if timepoints is None:
            timepoints = range(nd2.sizes["t"])
    else:
        if timepoints is None:
            timepoints = it.cycle(range(nd2.sizes["t"]))
    ts_iter = list(zip(timepoints, scale))
    if progress_bar is not None:
        ts_iter = progress_bar(ts_iter)
    animation = [
        mosaic_frame(
            get_frame_func,
            positions,
            input_dims,
            timepoint=t,
            scale=s,
            center=center,
            offset=offset,
            delayed=delayed,
        )
        for t, s in ts_iter
    ]
    return animation


def get_intensity_extrema(nd2, channels, v=0, step=10):
    extrema = {}
    for channel in channels:
        min_value = -1
        max_value = -1
        for t in range(0, nd2.sizes["t"], step):
            img = nd2.get_frame_2D(v=v, t=t, c=nd2.metadata["channels"].index(channel))
            if min_value == -1:
                min_value = img.min()
                # max_value = img.max()
                max_value = np.percentile(img, 99.9)
            else:
                min_value = min(min_value, img.min())
                # max_value = max(max_value, img.max())
                max_value = max(max_value, np.percentile(img, 99.9))
        extrema[channel] = (min_value, max_value)
    return extrema


def get_scaling_funcs(extrema):
    scaling_funcs = {}
    for channel, (min_value, max_value) in extrema.items():
        # careful! there's an unfortunate late-binding issue
        # SEE: https://stackoverflow.com/questions/1107210/python-create-function-in-a-loop-capturing-the-loop-variable
        # TODO: this should be a single clip...
        scaling_funcs[
            channel
        ] = lambda x, min_value=min_value, max_value=max_value: np.clip(
            (np.clip(x, min_value, max_value) - min_value) / (max_value - min_value),
            0,
            1,
        )
    return scaling_funcs


def export_video(ary, filename, fps=30, codec="h264", crf=22, tune="stillimage"):
    with av.open(filename, mode="w") as container:
        stream = container.add_stream(
            codec, rate=fps, options={"crf": str(crf), "tune": tune}
        )
        stream.width = ary[0].shape[1]
        stream.height = ary[0].shape[0]
        stream.pix_fmt = "yuv420p"
        for idx in range(len(ary)):
            img = np.round(255 * ary[idx]).astype(np.uint8)
            img = np.clip(img, 0, 255)
            frame = av.VideoFrame.from_ndarray(img, format="rgb24")
            for packet in stream.encode(frame):
                container.mux(packet)
        for packet in stream.encode():
            container.mux(packet)
