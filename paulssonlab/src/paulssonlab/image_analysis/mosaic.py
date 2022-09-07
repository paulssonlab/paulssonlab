import numpy as np
from matplotlib.colors import hex2color
import dask
import dask.array as da
import scipy.ndimage as ndi
import skimage
from skimage.transform import SimilarityTransform, warp
from PIL import Image, ImageDraw, ImageFont
import nd2reader
import av
import numba
import itertools as it
from functools import cache
from numbers import Number
from cytoolz import partial
from tqdm.auto import tqdm
from paulssonlab.image_analysis.workflow import (
    parse_nd2_metadata,
    get_position_metadata,
    get_filename_image_limits,
    get_nd2_frame,
)
from paulssonlab.image_analysis.blur import scipy_box_blur

# TODO: move to paulssonlab.util
def get_delayed(delayed):
    if delayed is True:
        return dask.delayed(pure=True)
    elif delayed is False:
        return lambda func, **kwargs: func
    else:
        return delayed


def colorize(imgs, hexcolors, scale=True):
    colors = [hex2color(hexcolor) for hexcolor in hexcolors]
    return _colorize(imgs, colors, scale=scale)


# @numba.njit
def _colorize(channel_imgs, colors, scale=True):
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
        scaled_img[:, :, np.newaxis] * np.array(color)
        for scaled_img, color in zip(scaled_imgs, colors)
    ]
    if not len(imgs_to_combine):
        return np.ones(channel_imgs[0].shape)  # white placeholder
    img = imgs_to_combine[0]
    for img2 in imgs_to_combine[1:]:
        img = 1 - (1 - img) * (1 - img2)
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


@cache
def fit_output_transform(input_width, input_height, output_width, output_height):
    width_ratio = input_width / output_width
    height_ratio = input_height / output_height
    scale = max(width_ratio, height_ratio)
    x = -(output_width - input_width / scale) / 2
    y = -(output_height - input_height / scale) / 2
    return SimilarityTransform(translation=(x, y)) + SimilarityTransform(scale=scale)


@cache
def fit_output_and_scale_transform(
    input_width,
    input_height,
    output_width,
    output_height,
    center_x,
    center_y,
    output_corner_x,
    output_corner_y,
    scale,
):
    transform = (
        fit_output_transform(input_width, input_height, output_width, output_height)
        + scale_around_center(1 / scale, (input_width / 2, input_height / 2))
        + SimilarityTransform(
            translation=(
                center_x - input_width / 2,
                center_y - input_height / 2,
            )
        )
        + SimilarityTransform(translation=(output_corner_x, output_corner_y))
    )
    input_ul = transform.inverse((0, 0))[0]
    # TODO: off-by-one?
    input_lr = transform.inverse((input_width - 1, input_height - 1))[0]
    output_ul = (0, 0)
    # TODO: off-by-one?
    output_lr = (output_width - 1, output_height - 1)
    visible = rectangles_intersect(input_ul, input_lr, output_ul, output_lr)
    return transform, visible


def mosaic_frame(
    get_frame_func,
    channels,
    colors,
    positions,
    input_dims,
    *,
    timepoint=None,
    center=None,
    offset=None,
    scale=1,
    output_dims=(1024, 1024),
    anti_aliasing=True,
    max_gaussian_sigma=4,
    scaling_funcs=None,
    delayed=False,
    dtype=np.float32,
):
    delayed = get_delayed(delayed)
    if len(channels) != len(colors):
        raise ValueError("number of channels and colors must be equal")
    if center is not None and offset is not None:
        raise ValueError("can only specify at most one of center, offset")
    if center is None:
        columns = positions["x_idx"].max() - positions["x_idx"].min() + 1
        rows = positions["y_idx"].max() - positions["y_idx"].min() + 1
        center = np.array([input_dims[0] * columns / 2, input_dims[1] * rows / 2])
    if offset is not None:
        center += np.array(offset)
    all_channel_imgs = [[] for _ in range(len(channels))]
    for (filename, pos_num), position in positions.iterrows():
        frame_transform, visible = fit_output_and_scale_transform(
            *input_dims,
            *output_dims,
            *center,
            -input_dims[0] * position["x_idx"],
            -input_dims[1] * position["y_idx"],
            scale,
        )
        if visible:
            for channel, channel_imgs in zip(channels, all_channel_imgs):
                img = delayed(get_frame_func)(pos_num, channel, timepoint)
                if scaling_funcs:
                    img = delayed(scaling_funcs[channel])(img)
                if anti_aliasing:
                    # taken from skimage.transform.resize, with an added factor of 1/4
                    # (because original formula gave sigmas that were too large)
                    # this was not based on any well-founded heuristic
                    # SEE: https://github.com/scikit-image/scikit-image/pull/2802
                    anti_aliasing_sigma = max(0, (frame_transform.scale - 1) / 8)
                    if anti_aliasing_sigma > 0:
                        if anti_aliasing_sigma <= max_gaussian_sigma:
                            blur_func = ndi.gaussian_filter
                        else:
                            # TODO: use stack blur instead?
                            blur_func = scipy_box_blur
                        img = delayed(blur_func)(
                            img,
                            anti_aliasing_sigma,
                            mode="nearest",
                        )
                img = delayed(warp)(
                    img, frame_transform, output_shape=output_dims[::-1], order=1
                )
                img = da.from_delayed(img, output_dims[::-1], dtype=dtype)
                channel_imgs.append(img)
    channel_composite_imgs = [
        da.stack(channel_imgs).sum(axis=0) for channel_imgs in all_channel_imgs
    ]
    # TODO: implement delayed=False without dask.array? or just call .compute?
    output_img = delayed(colorize)(
        channel_composite_imgs, colors, scale=(not scaling_funcs)
    )
    return output_img


def _composite_rgba(rgb, rgba):
    alpha = rgba[:, :, -1, np.newaxis]
    return alpha * rgb + (1 - alpha) * rgba[:, :, :-1]


def square_overlay(
    frame,
    timepoint,
    scale,
    unit_scale=80,
    unit=1,
    unit_width=0.5,
    unit_noun=("cell", "cells"),
    factor=10,
    caption_position="inset",
    font=None,
    font_size=60,
    text_padding=20,
    line_width=3,
    color=(1, 1, 1),
):
    color = (*np.array(color) * 255, 0)
    img = Image.new("RGBA", frame.shape[:-1], (0, 0, 0, 255))
    draw = ImageDraw.Draw(img)
    min_dim = min(*frame.shape[:-1])
    V = (min_dim / scale) ** 2
    U = (min_dim * unit_width) ** 2 / (unit * unit_scale**2)
    n = np.floor((np.log(V) - np.log(U)) / np.log(factor))
    count = factor**n
    width = scale / min_dim * np.sqrt(U * count)
    half_width_px = min_dim * width / 2
    center = np.array([frame.shape[1] / 2, frame.shape[0] / 2])
    delta = np.array([half_width_px, -half_width_px])
    draw.rectangle(
        [tuple(center - delta), tuple(center + delta)], outline=color, width=line_width
    )
    match caption_position:
        case "inset":
            if count == 1:
                noun = unit_noun[0]
            else:
                noun = unit_noun[1]
            caption = f"{int(count)} {noun}"
            scaled_font_size = int(np.ceil(font_size * width))
            text_padding = np.array([text_padding, -scaled_font_size - text_padding])
            draw.text(
                tuple(center - delta + text_padding),
                caption,
                font=font.font_variant(size=scaled_font_size),
                fill=color,
            )
        case _:
            raise ValueError("caption_position not recognized")
    return _composite_rgba(frame, np.asarray(img) / 255)


def mosaic_animate_scale(
    filename,
    scale=1,
    timepoints=None,
    center=None,
    offset=None,
    width=1024,
    height=1024,
    channels=None,
    channel_to_color=None,
    scaling_funcs=None,
    overlay_func=None,
    delayed=True,
    progress_bar=tqdm,
):
    if channels is None:
        raise ValueError("must specify channels")
    if channel_to_color is None:
        raise ValueError("must specify channel_to_color")
    if overlay_func:

        def frame_func(*args, timepoint=None, scale=None, **kwargs):
            frame = mosaic_frame(*args, timepoint=timepoint, scale=scale, **kwargs)
            return delayed(overlay_func)(frame, timepoint, scale)

    else:
        frame_func = mosaic_frame
    delayed = get_delayed(delayed)
    nd2 = nd2reader.ND2Reader(filename)
    nd2s = {filename: nd2 for filename in (filename,)}
    metadata = {
        nd2_filename: parse_nd2_metadata(nd2) for nd2_filename, nd2 in nd2s.items()
    }
    positions = get_position_metadata(metadata)
    image_limits = get_filename_image_limits(metadata)
    get_frame_func = partial(
        get_nd2_frame,
        filename,
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
    colors = [channel_to_color[channel] for channel in channels]
    ts_iter = list(zip(timepoints, scale))
    if progress_bar is not None:
        ts_iter = progress_bar(ts_iter)
    animation = [
        frame_func(
            get_frame_func,
            channels,
            colors,
            positions,
            input_dims,
            timepoint=t,
            scale=s,
            center=center,
            offset=offset,
            scaling_funcs=scaling_funcs,
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
                max_value = np.percentile(img, 99.9)
            else:
                min_value = min(min_value, img.min())
                max_value = max(max_value, np.percentile(img, 99.9))
        extrema[channel] = (min_value, max_value)
    return extrema


def get_scaling_funcs(extrema):
    scaling_funcs = {}
    for channel, (min_value, max_value) in extrema.items():
        # careful! there's an unfortunate late-binding issue
        # SEE: https://stackoverflow.com/questions/1107210/python-create-function-in-a-loop-capturing-the-loop-variable
        scaling_funcs[channel] = lambda x, min_value=min_value, max_value=max_value: (
            np.clip(x, min_value, max_value) - min_value
        ) / (max_value - min_value)
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
