import holoviews as hv
import imageio.v3 as iio
import IPython.display
import numpy as np


def display_image(img, scale=False, downsample=1, format="jpg"):
    if downsample != 1:
        img = img[::downsample, ::downsample, ...]
    if scale is True:
        img_min = np.nanmin(img)
        img = (img - img_min) / (np.nanmax(img) - img_min)
    elif scale:
        img = img - np.nanmin(img)
        img = img / np.nanpercentile(img, scale * 100)
    img = (np.clip(img, 0, 1) * 255).astype(np.uint8)
    if format.lower() in ("jpg", "jpeg"):
        bytes_ = iio.imwrite("<bytes>", img, extension=".jpeg", quality=95)
    elif format.lower() == "png":
        bytes_ = iio.imwrite("<bytes>", img, extension=".png")
    else:
        raise ValueError(f"unknown format: {format}")
    return IPython.display.Image(bytes_)


def RevImage(img, **kwargs):
    return _RevImage(hv.Image, img, **kwargs)


def _RevImage(cls, img, scale=False, **kwargs):
    plot = cls(
        img[::-1],
        **{
            "bounds": (0, 0, img.shape[1], img.shape[0]),
            **kwargs,
        },
    )
    if scale is True:
        img_min = np.nanmin(img)
        img_max = np.nanmax(img)
        plot = plot.redim.range(z=(img_min, img_max))
    elif scale:
        img_min = np.nanmin(img)
        img_max = np.nanpercentile(img, scale * 100)
        plot = plot.redim.range(z=(img_min, img_max))
    return plot.opts(aspect=img.shape[1] / img.shape[0], invert_yaxis=True)


def RevRGB(img, **kwargs):
    return _RevImage(hv.RGB, img, **kwargs)
