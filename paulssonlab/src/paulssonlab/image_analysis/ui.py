import holoviews as hv
import numpy as np


def RevImage(img, **kwargs):
    return _RevImage(hv.Image, img, **kwargs)


def _RevImage(cls, img, **kwargs):
    # TODO: bug in holoviews/bokeh aspect ratio handling, you need data_aspect=(height/width)**2 instead of 1
    # TODO: I think datashader is doing this correctly, but the original Image display screws up the frame dimensions
    # return cls(img[::-1], bounds=(0,0,img.shape[1],img.shape[0])).options(invert_yaxis=True, data_aspect=(img.shape[1]/img.shape[0])**2, width=400)#, responsive=False)
    return cls(img[::-1], bounds=(0, 0, img.shape[1], img.shape[0])).options(
        aspect=img.shape[1] / img.shape[0], frame_width=250
    )  # , responsive=False)


def RevRGB(img, **kwargs):
    return _RevImage(hv.RGB, img, **kwargs)
