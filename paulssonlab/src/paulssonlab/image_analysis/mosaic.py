import numpy as np
from matplotlib.colors import hex2color
import dask
import dask.array as da
import distributed
import scipy.ndimage as ndi
import skimage
from skimage.transform import SimilarityTransform, warp
from PIL import Image, ImageDraw, ImageFont
import cv2
import nd2reader
import av
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
    input_scale = fixed_aspect_scale(
        *input_dims, output_dims[0] * scale, output_dims[1] * scale
    )
    rescaled_input_dims = np.ceil(np.array(input_dims) * input_scale).astype(np.int_)
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
                img = delayed(get_frame_func)(pos_num, channel, timepoint)
                if scaling_funcs:
                    img = delayed(scaling_funcs[channel])(img)
                img = delayed(cv2.resize)(
                    img, rescaled_input_dims, interpolation=cv2.INTER_AREA
                )
                # TODO: cv2.INTER_AREA is not implemented for cv2.warpAffine,
                # so we resize first, then translate/rotate
                img = delayed(cv2.warpAffine)(
                    img,
                    frame_transform.params[:2, :],
                    output_dims,
                    # flags=cv2.INTER_AREA + cv2.WARP_INVERSE_MAP,
                    flags=(cv2.INTER_LANCZOS4 + cv2.WARP_INVERSE_MAP),
                )
                img = delayed(np.clip)(img, 0, 1)  # LANCZOS4 outputs values beyond 0..1
                img = da.from_delayed(img, output_dims[::-1], dtype=dtype)
                channel_imgs.append(img)
    if not all_channel_imgs:
        raise ValueError("no positions visible")
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
    return (1 - alpha) * rgb + alpha * rgba[:, :, :-1]


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
    noun=("cell", "cells"),
    factor=10,
    font=None,
    font_size=60,
    text_padding=20,
    line_width=3,
    color=(1, 1, 1),
    alpha_width=60,
    alpha_transition=10,
    text_alpha_transition=5,
):
    img = Image.new("RGBA", frame.shape[:-1:][::-1], (0, 0, 0, 255))
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
    for n in range(min_n, viewport_n + 2):
        count = factor**n
        geometry_scale = scale / min_dim * np.sqrt(min_U * fit_factor ** (n - min_n))
        half_width_px = min_dim * geometry_scale / 2
        alpha = 1 / (1 + np.exp(-(half_width_px - alpha_width) / alpha_transition))
        if alpha < 0.05:
            continue
        rgba_color = np.array((*np.array(color) * 255, np.ceil(alpha * 255))).astype(
            np.uint8
        )
        center = np.array([frame.shape[1] / 2, frame.shape[0] / 2])
        delta = np.array([half_width_px, half_width_px])
        draw.rectangle(
            [tuple(center - delta), tuple(center + delta)],
            outline=tuple(rgba_color),
            width=line_width,
            # width=int(np.ceil(line_width * geometry_scale)),
        )
        caption = f"{int(count)} {noun[0] if count == 1 else noun[1]}"
        scaled_font_size = int(np.ceil(font_size * geometry_scale))
        text_padding_ary = np.array([text_padding, -scaled_font_size - text_padding])
        # if width >= alpha_width:
        #     text_alpha = alpha
        # else:
        #     text_alpha = 1 - 1 / (
        #         1 + np.exp(-(half_width_px - alpha_width) / text_alpha_transition)
        #     )
        # if text_alpha > 0.95:
        #     continue
        # text_rgba_color = (*np.array(color) * 255, int(np.ceil(text_alpha * 255)))
        text_rgba_color = rgba_color
        # print("!", scaled_font_size, tuple(center - delta + text_padding_ary))
        if scaled_font_size > 3:
            text_position = (
                center + np.array([-half_width_px, half_width_px]) + text_padding_ary
            )
            # draw.text(
            #     tuple(text_position),
            #     caption,
            #     font=font.font_variant(size=scaled_font_size),
            #     fill=tuple(text_rgba_color),
            # )
            #####
            text_img = Image.new("RGBA", frame.shape[:-1:][::-1], (0, 0, 0, 0))
            text_draw = ImageDraw.Draw(text_img)
            text_padding_ary2 = np.array([text_padding, -font_size - text_padding])
            text_position2 = (
                center + np.array([-min_dim / 2, min_dim / 2]) + text_padding_ary2
            )
            text_rgba_color2 = np.array(
                (*np.array(color) * 255, int(np.ceil(alpha * 255)))
            )
            # text_rgba_color2 = np.array((*np.array(color) * 255, 255))
            text_draw.text(
                tuple(text_position2),
                caption,
                font=font.font_variant(size=font_size),
                fill=tuple(text_rgba_color2),
            )
            text_img_dims = np.ceil(np.array(frame.shape[:-1]) * geometry_scale).astype(
                np.int32
            )
            text_img_resized = cv2.resize(
                np.asarray(text_img) / 255,
                text_img_dims,
                interpolation=cv2.INTER_AREA,
            )
            # transform = SimilarityTransform(translation=np.array(frame.shape[:-1]) / 2)
            transform = SimilarityTransform(
                translation=(np.array(frame.shape[:-1]) - text_img_dims) // 2
            )
            text_img_translated = cv2.warpAffine(
                text_img_resized,
                transform.params[:2, :],
                frame.shape[:-1],
                flags=(cv2.INTER_LANCZOS4),
            )
            # transform = SimilarityTransform(scale=width)
            # img2_rescaled = warp(np.asarray(img2) / 255, transform)
            # print(img2_rescaled.dtype, img2_rescaled.shape)
            img.alpha_composite(
                Image.fromarray((text_img_translated * 255).astype(np.uint8))
            )
            # img.alpha_composite(img2)
    # return np.asarray(img)[:, :, :3] / 255
    # return _composite_rgba(frame, np.asarray(text_ary) / 255)
    return _composite_rgba(frame, np.asarray(img) / 255)


def mosaic_animate_scale(
    filename,
    scale=1,
    timepoints=None,
    center=None,
    offset=None,
    rotation=None,
    output_dims=(1024, 1024),
    channels=None,
    channel_to_color=None,
    scaling_funcs=None,
    overlay_func=None,
    overlay_only=False,
    positions_func=None,
    delayed=True,
    progress_bar=tqdm,
):
    if channels is None:
        raise ValueError("must specify channels")
    if channel_to_color is None:
        raise ValueError("must specify channel_to_color")
    if overlay_func:

        def frame_func(*args, timepoint=None, scale=None, **kwargs):
            if overlay_only:
                frame = np.zeros((*output_dims, 3), dtype=np.float32)
            else:
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
    if positions_func is not None:
        positions = positions_func(positions)
    if callable(rotation):
        rotation = rotation(positions)
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
            rotation=rotation,
            scaling_funcs=scaling_funcs,
            output_dims=output_dims,
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


def export_video(
    ary, filename, fps=30, codec="h264", crf=22, tune="stillimage", progress_bar=tqdm
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
            if isinstance(img, distributed.Future):
                img = img.result()
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
