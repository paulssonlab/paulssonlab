import numpy as np
from scipy.special import binom
import freetype
import gdspy as g

DEFAULT_FACE = freetype.Face("Vera.ttf")


def bernstein(n, i, t):
    return binom(n, i) * (t**i) * ((1 - t) ** (n - i))


def bezier(points, t):
    # CF: https://github.com/vicgalle/machine-learning/blob/master/bezier.py
    n = len(points) - 1
    return (
        bernstein(n, np.arange(n + 1)[:, np.newaxis], t)[..., np.newaxis]
        * points[:, np.newaxis]
    ).sum(axis=0)


def Text2(
    text,
    size,
    position=(0, 0),
    horizontal=True,
    angle=0,
    x_reflection=False,
    layer=0,
    datatype=0,
    face=DEFAULT_FACE,
):
    face.set_char_size(48 * 64)
    face.load_char("S")
    slot = face.glyph
    outline = slot.outline
