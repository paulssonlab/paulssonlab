import itertools as it
from collections.abc import Sequence
from functools import cache, partial
from numbers import Number

import av
import cv2
import dask
import dask.array as da
import numpy as np
import scipy.ndimage as ndi
import skimage
import zarr
from dask.distributed import Future
from PIL import Image, ImageDraw, ImageFont
from skimage.transform import SimilarityTransform, warp
from tqdm.auto import tqdm

from paulssonlab.image_analysis.blur import scipy_box_blur
from paulssonlab.image_analysis.image import colorize
from paulssonlab.image_analysis.util import get_delayed
from paulssonlab.image_analysis.workflow import (
    get_filename_image_limits,
    get_nd2_frame,
    get_nd2_reader,
    get_position_metadata,
)
from paulssonlab.io.metadata import parse_nd2_metadata


def composite_rgb_rgba(rgb, rgba):
    alpha = rgba[:, :, -1, np.newaxis]
    return (1 - alpha) * rgb + alpha * rgba[:, :, :-1]


# TODO: unused
def composite_mean(stack, fillval=None):
    mean = da.average(stack, axis=0, weights=(stack != 0))
    if fillval is not None:
        mean = da.where(da.isnan(mean), fillval, mean)
    return mean


# TODO: unused
# use via composite_first(da.stack(channel_imgs).rechunk((-1, "auto", "auto")))
def composite_first(stack):
    # TODO: this doesn't work because second argument needs to be a numpy array, not delayed dask array
    # return da.choose(stack, (stack != 0).argmax(axis=0))
    return da.map_blocks(
        lambda ary, idxs: np.take_along_axis(ary, idxs[np.newaxis, ...], axis=0)[0],
        stack,
        (stack != 0).argmax(axis=0),
        dtype=stack.dtype,
        name=f"composite_first-{tokenize(stack)}",
        drop_axis=0,
    )


def rectangles_intersect(ul1, lr1, ul2, lr2):
    return not (
        (ul1[0] > lr2[0]) or (lr1[0] < ul2[0]) or (ul1[1] > lr2[1]) or (lr1[1] < ul2[1])
    )


def scale_around_center(scale, center):
    x, y = center
    return (
        SimilarityTransform(translation=(-x, -y))
        + SimilarityTransform(scale=scale)
        + SimilarityTransform(translation=(scale * x, scale * y))
    )


def rotate_around_center(rotation, width, height):
    # FROM: https://stackoverflow.com/a/36045144
    center = (np.array([width, height]) - 1) / 2
    return (
        SimilarityTransform(translation=-center)
        + SimilarityTransform(rotation=rotation)
        + SimilarityTransform(translation=center)
    )


def fixed_aspect_scale(input_width, input_height, output_width, output_height):
    width_ratio = output_width / input_width
    height_ratio = output_height / input_height
    scale = min(width_ratio, height_ratio)
    return scale


@cache
def transform_to_viewport(
    input_width,
    input_height,
    output_width,
    output_height,
    center_x,
    center_y,
    output_corner_x,
    output_corner_y,
    rotation=None,
):
    translation = SimilarityTransform(
        translation=(
            center_x - output_width / 2 + output_corner_x,
            center_y - output_height / 2 + output_corner_y,
        ),
    )
    rotation = rotate_around_center(rotation, output_width, output_height)
    transform = translation + rotation
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
    rotation=None,
    scale=1,
    output_dims=(1024, 1024),
    scaling_funcs=None,
    dark=None,
    flats=None,
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
    input_scale = scale * fixed_aspect_scale(
        *input_dims, output_dims[0], output_dims[1]
    )
    rescaled_input_dims = np.ceil(np.array(input_dims) * input_scale).astype(np.uint32)
    for (filename, pos_num), position in positions.iterrows():
        frame_transform, visible = transform_to_viewport(
            *rescaled_input_dims,
            *output_dims,
            *(center * input_scale),
            -input_dims[0] * position["x_idx"] * input_scale,
            -input_dims[1] * position["y_idx"] * input_scale,
            rotation=rotation,
        )
        if visible:
            for channel, channel_imgs in zip(channels, all_channel_imgs):
                img = delayed(get_frame_func)(
                    pos_num,
                    channel,
                    timepoint,
                    dark=dark,
                    flat=(flats or {}).get(channel),
                )
                if scaling_funcs:
                    img = delayed(scaling_funcs[channel])(img)
                if scale < 1:
                    # TODO: cv2.INTER_AREA is not implemented for cv2.warpAffine,
                    # so for downscaling (scale < 1) we resize first, then translate/rotate
                    img = delayed(cv2.resize)(
                        img, rescaled_input_dims, interpolation=cv2.INTER_AREA
                    )
                    frame_transform_with_scaling = frame_transform
                else:
                    scale_transform = scale_around_center(
                        1 / input_scale, np.array(input_dims) / 2 - 1
                    )
                    frame_transform_with_scaling = frame_transform + scale_transform
                img = delayed(cv2.warpAffine)(
                    img,
                    frame_transform_with_scaling.params[:2, :],
                    output_dims,
                    flags=(cv2.INTER_LANCZOS4 + cv2.WARP_INVERSE_MAP),
                )
                img = delayed(np.clip)(img, 0, 1)  # LANCZOS4 outputs values beyond 0..1
                # TODO: make this support the delayed=False case
                img = da.from_delayed(img, output_dims[::-1], dtype=dtype)
                channel_imgs.append(img)
    if not all_channel_imgs:
        raise ValueError("no positions visible")
    channel_composite_imgs = [
        da.stack(channel_imgs).max(axis=0) for channel_imgs in all_channel_imgs
    ]
    output_img = delayed(colorize)(
        channel_composite_imgs, colors, scale=(not scaling_funcs)
    )
    return output_img


def square_overlay(
    frame,
    timepoint,
    scale,
    min_scale=80,
    min_n=0,
    min_width=0.5,
    max_scale=0.1,
    max_n=5,
    max_width=0.9,
    n_range=None,
    noun=("cell", "cells"),
    factor=10,
    font=None,
    font_size=180,
    text_padding=40,
    line_width=20,
    color=(1, 1, 1),
    alpha_width=60,
    alpha_transition=10,
    text_alpha_transition=5,
):
    if n_range is None:
        n_range = (min_n, max_n)
    output_dims = np.array(frame.shape[:-1:][::-1])
    img = Image.new("RGBA", tuple(output_dims), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)
    min_dim = min(*frame.shape[:-1])
    V = (min_dim / scale) ** 2
    # min_U, max_U are the number of pixels
    min_U = (min_dim * min_width) ** 2 / min_scale**2
    max_U = (min_dim * max_width) ** 2 / max_scale**2
    # calculate fit_factor so that squares at min_scale and max_scale
    # result in counts of factor**min_n, factor**max_n, respectively
    fit_factor = (max_U / min_U) ** (1 / (max_n - min_n))
    viewport_n = int(np.floor((np.log(V) - np.log(min_U)) / np.log(fit_factor)))
    for n in range(n_range[0], min(viewport_n + 2, n_range[1])):
        count = factor**n
        geometry_scale = scale / min_dim * np.sqrt(min_U * fit_factor ** (n - min_n))
        half_width_px = min_dim * geometry_scale / 2
        alpha = 1 / (1 + np.exp(-(half_width_px - alpha_width) / alpha_transition))
        if alpha < 0.05:
            continue
        rgba_color = np.array((*np.array(color) * 255, np.ceil(alpha * 255))).astype(
            np.uint8
        )
        center = output_dims / 2
        delta = np.array([half_width_px, half_width_px])
        draw.rectangle(
            [tuple(center - delta), tuple(center + delta)],
            outline=tuple(rgba_color),
            # width=line_width,
            width=int(np.ceil(line_width * geometry_scale)),
        )
        caption = f"{int(count)} {noun[0] if count == 1 else noun[1]}"
        scaled_font_size = int(np.ceil(font_size * geometry_scale))
        text_rgba_color = rgba_color
        if scaled_font_size > 3:
            text_img = Image.new("RGBA", tuple(output_dims), (0, 0, 0, 0))
            text_draw = ImageDraw.Draw(text_img)
            text_padding_ary = np.array([text_padding, -font_size - text_padding])
            text_position = (
                center + np.array([-min_dim / 2, min_dim / 2]) + text_padding_ary
            )
            text_draw.text(
                tuple(text_position),
                caption,
                font=font.font_variant(size=font_size),
                fill=tuple(rgba_color),
            )
            text_img_dims = np.ceil(np.array(output_dims) * geometry_scale).astype(
                np.int32
            )
            text_img_resized = cv2.resize(
                np.asarray(text_img) / 255,
                text_img_dims,
                interpolation=cv2.INTER_AREA,
            )
            del text_draw
            del text_img
            transform = SimilarityTransform(
                translation=(np.array(output_dims) - text_img_dims) // 2
            )
            text_img_translated = cv2.warpAffine(
                text_img_resized,
                transform.params[:2, :],
                output_dims,
                flags=cv2.INTER_LINEAR,
            )
            del text_img_resized
            img.alpha_composite(
                Image.fromarray((text_img_translated * 255).astype(np.uint8))
            )
            del text_img_translated
    return composite_rgb_rgba(frame, np.asarray(img) / 255)


def get_nd2_metadata(filename):
    nd2 = get_nd2_reader(filename)
    nd2s = {filename: nd2 for filename in (filename,)}
    metadata = {
        nd2_filename: parse_nd2_metadata(nd2) for nd2_filename, nd2 in nd2s.items()
    }
    positions = get_position_metadata(metadata)
    image_limits = get_filename_image_limits(metadata)
    input_dims = (
        image_limits[filename][0][1] + 1,
        image_limits[filename][1][1] + 1,
    )
    return positions, input_dims


def mosaic_animate_scale(
    get_frame_func,
    scale=1,
    timepoints=None,
    positions=None,
    center=None,
    offset=None,
    rotation=None,
    input_dims=None,
    output_dims=(1024, 1024),
    channels=None,
    channel_to_color=None,
    scaling_funcs=None,
    overlay_func=None,
    overlay_only=False,
    dark=None,
    flats=None,
    delayed=True,
    progress_bar=tqdm,
):
    delayed = get_delayed(delayed)
    if positions is None:
        raise ValueError("must specify positions")
    if input_dims is None:
        raise ValueError("must specify input_dims")
    if channels is None:
        raise ValueError("must specify channels")
    if channel_to_color is None:
        raise ValueError("must specify channel_to_color")
    if overlay_func:

        def frame_func(*args, timepoint=None, scale=None, **kwargs):
            if overlay_only:
                frame = delayed(np.zeros)((*output_dims[::-1], 3), dtype=np.float32)
            else:
                frame = mosaic_frame(*args, timepoint=timepoint, scale=scale, **kwargs)
            return delayed(overlay_func)(frame, timepoint, scale)

    else:
        frame_func = mosaic_frame
    if isinstance(channels, Sequence) and isinstance(channels[0], str):
        channels = it.repeat(channels)
    tsc_iter = list(zip(timepoints, scale, channels))
    if progress_bar is not None:
        tsc_iter = progress_bar(tsc_iter)
    animation = [
        frame_func(
            get_frame_func,
            frame_channels,
            [channel_to_color[channel] for channel in frame_channels],
            positions,
            input_dims,
            timepoint=t,
            scale=s,
            center=center,
            offset=offset,
            rotation=rotation,
            scaling_funcs=scaling_funcs,
            output_dims=output_dims,
            dark=dark,
            flats=flats,
            delayed=delayed,
        )
        for t, s, frame_channels in tsc_iter
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


def max_scale(positions, input_dims, output_dims, offset=None, mode=None):
    input_dims = np.asarray(input_dims)
    output_dims = np.asarray(output_dims)
    if offset is None:
        offset = np.array([0, 0])
    else:
        offset = np.asarray(offset)
    columns = positions["x_idx"].max() - positions["x_idx"].min() + 1
    rows = positions["y_idx"].max() - positions["y_idx"].min() + 1
    input_scale = fixed_aspect_scale(*input_dims, output_dims[0], output_dims[1])
    scale = (
        output_dims
        / (
            np.array([columns, rows]) * input_dims
            + 2 * np.array([-1, 1])[:, np.newaxis] * offset
        )
        / input_scale
    )
    if mode is None:
        return scale
    elif mode == "fit":
        return scale.min()
    elif mode == "fill":
        return scale.max()
    elif mode == "fit_horizontal":
        return scale[:, 0].min()
    elif mode == "fill_horizontal":
        return scale[:, 0].max()
    elif mode == "fit_vertical":
        return scale[:, 1].min()
    elif mode == "fill_vertical":
        return scale[:, 1].max()
    else:
        raise ValueError(f"unknown mode: {mode}")


def export_video(
    ary,
    filename,
    fps=30,
    downsample=None,
    codec="h264",
    crf=22,
    tune="stillimage",
    progress_bar=tqdm,
):
    with av.open(str(filename), mode="w") as container:
        stream = container.add_stream(
            codec, rate=fps, options={"crf": str(crf), "tune": tune}
        )
        initialized = False
        idxs = range(len(ary))
        if progress_bar is not None:
            idxs = progress_bar(idxs)
        for idx in idxs:
            img = ary[idx]
            if isinstance(img, Future):
                img = img.result()
            if downsample is not None:
                img = img[::downsample, ::downsample, ...]
            if not initialized:
                stream.width = img.shape[1]
                stream.height = img.shape[0]
                stream.pix_fmt = "yuv420p"
                initialized = True
            img = np.round(255 * img).astype(np.uint8)
            img = np.clip(img, 0, 255)
            frame = av.VideoFrame.from_ndarray(img, format="rgb24")
            for packet in stream.encode(frame):
                container.mux(packet)
        for packet in stream.encode():
            container.mux(packet)


def write_to_zarr(filename, frame, frame_num, num_frames, dtype=np.float32):
    arr = zarr.open_array(
        filename,
        mode="a",
        shape=(num_frames, *frame.shape),
        chunks=(1, *frame.shape),
        dtype=dtype,
    )
    arr[frame_num, ...] = frame
