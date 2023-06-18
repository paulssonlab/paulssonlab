import holoviews as hv
import numpy as np


def RevImage(img, **kwargs):
    return _RevImage(hv.Image, img, **kwargs)


def _RevImage(cls, img, scale=False, **kwargs):
    plot = cls(
        img[::-1, ...],
        **{
            "bounds": (0, img.shape[0], img.shape[1], 0),
            "extents": (0, img.shape[0], img.shape[1], 0),
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


def overlay_inverted_yaxis(overlay):
    return overlay.opts(
        hv.opts.Overlay(invert_yaxis=True), hv.opts.Image(invert_yaxis=False)
    )
