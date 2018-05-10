import numpy as np
from scipy.special import binom
import freetype
import gdspy as g
from geometry import Cell
from util import memoize

DEFAULT_FACE = freetype.Face("Vera.ttf")
POINTS_PER_SEGMENT = 10


def bernstein(n, i, t):
    return binom(n, i) * (t**i) * ((1 - t) ** (n - i))


def bezier(points, t):
    # SEE: https://github.com/vicgalle/machine-learning/blob/master/bezier.py
    n = len(points) - 1
    return (
        bernstein(n, np.arange(n + 1)[:, np.newaxis], t)[..., np.newaxis]
        * points[:, np.newaxis]
    ).sum(axis=0)


def render_character(
    char, size, face=DEFAULT_FACE, points_per_segment=POINTS_PER_SEGMENT
):
    # SEE: https://www.freetype.org/freetype2/docs/glyphs/glyphs-6.html
    face.set_char_size(size)
    face.load_char(char)
    outline = face.glyph.outline
    contours = np.array(outline.contours)
    polygons = []
    for start, end in zip(np.concatenate(((0,), contours[:-1] + 1)), contours + 1):
        points = outline.points[start:end]
        tags = outline.tags[start:end]
        segments = []
        for j in range(len(points)):
            if (
                tags[j] & freetype.FT_CURVE_TAG_CONIC
                and j > 1
                and tags[j - 1] & freetype.FT_CURVE_TAG_CONIC
            ):
                # create virtual on point at midpoint of points[j-1] and points[j]
                midpoint = (points[j - 1] + points[j]) / 2
                segments[-1].append(midpoint)
                segments.append([midpoint, points[j]])
            else:
                if j == 0 and tags[0] & freetype.FT_CURVE_TAG_CONIC:
                    # first point is off point, so segment should start at last point in contour
                    if tags[-1] & freetype.FT_CURVE_TAG_CONIC:
                        # last point is also off point, so segment should start and end with
                        # midpoint of first and last points
                        midpoint = (points[0] + points[-1]) / 2
                        segments.append([midpoint, points[0]])
                    else:
                        segments.append([points[-1], points[0]])
                if j != 0:
                    segments[-1].append(points[j])
                if tags[j] & freetype.FT_CURVE_TAG_ON:
                    segments.append([points[j]])
        if (
            tags[0] & freetype.FT_CURVE_TAG_CONIC
            and tags[-1] & freetype.FT_CURVE_TAG_CONIC
        ):
            # first and last points are off points, so last on point is midpoint of first and last points
            midpoint = (points[0] + points[-1]) / 2
            segments[-1].append[midpoint]
        else:
            # close last segment with first point
            segments[-1].append(points[0])
        polygon_segments = []
        # do not repeat last point, gdspy auto-closes polygons (linking last and first point)
        for segment in segments:
            if len(segment) == 2:
                polygon_segments.append((segment[0],))
            elif len(segment) == 3 or len(segment) == 4:
                polygon_segments.append(
                    bezier(np.array(segment), np.linspace(0, 1, points_per_segment))[
                        :-1
                    ]
                )
            else:
                raise ValueError(
                    "unexpected number of points per segment {}".format(len(segment))
                )
        polygon = np.concatenate(polygon_segments)
        polygons.append(polygon)
    metadata = {
        attr: getattr(face.glyph.metrics, attr)
        for attr in ("vertAdvance", "horiAdvance")
    }
    return polygons, metadata


@memoize
def character(
    char, size, face=DEFAULT_FACE, points_per_segment=POINTS_PER_SEGMENT, layer=0
):
    cell = Cell("~{}~s{}l{}".format(char, size, layer))
    polygons, metadata = render_character(
        char, size, face=face, points_per_segment=points_per_segment
    )
    cell.add(g.PolygonSet(polygons, layer=layer))
    return cell, metadata


def Text(
    text,
    size,
    position=(0, 0),
    horizontal=True,
    angle=0,
    layer=0,
    datatype=0,
    face=DEFAULT_FACE,
    points_per_segment=POINTS_PER_SEGMENT,
):
    if not horizontal or angle != 0:
        raise NotImplementedError
    position = np.array(position)
    offset = np.zeros(2)
    face.set_char_size(int(size * 64))
    linespace = face.size.height / 64
    face.load_char(" ")
    space_advance = face.glyph.metrics.horiAdvance / 64
    cells = []
    for char in text:
        if char == "\n":
            offset[0] = 0
            offset[1] -= linespace
            continue
        elif char == " ":
            face.load_char(char)
            offset[0] += space_advance
            continue
        char_cell, metadata = character(
            char, size, face=face, points_per_segment=points_per_segment
        )
        cells.append(
            g.CellReference(
                char_cell, position + offset, x_reflection=True, rotation=180
            )
        )
        offset[0] += metadata["horiAdvance"]
    return cells
