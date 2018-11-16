import numpy as np
from scipy.special import binom
import scipy.spatial
import freetype
import gdspy as g
from geometry import Cell, fast_boolean
from util import memoize, make_hashable
from functools import reduce, partial

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


def render_character(char, face=DEFAULT_FACE, points_per_segment=POINTS_PER_SEGMENT):
    # SEE: https://www.freetype.org/freetype2/docs/glyphs/glyphs-6.html
    face.load_char(char, freetype.FT_LOAD_NO_SCALE)
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
    metrics = {
        attr: getattr(face.glyph.metrics, attr)
        for attr in ("vertAdvance", "horiAdvance")
    }
    return polygons, metrics


def cut_out_hole(exterior, interior, layer=0, datatype=0):
    dists = scipy.spatial.distance.cdist(exterior, interior)
    idx_ext, idx_int = np.unravel_index(dists.argmin(), dists.shape)
    # exterior = np.concatenate((exterior[:idx_ext+1],
    #                               interior[:idx_int+1:-1],
    #                               interior[idx_int::-1],
    #                               exterior[idx_ext:]))
    result = np.concatenate(
        (
            exterior[: idx_ext + 1],
            interior[idx_int:],
            interior[: idx_int + 1],
            exterior[idx_ext:],
        )
    )
    return result


def cut_out_holes(polygons, layer=0, datatype=0):
    top_level_polygons = []
    used_polygons = set()
    for polygon_ext in polygons:
        polygons_to_cut = [polygon_ext]
        for polygon_int in polygons:
            if polygon_ext is polygon_int:
                continue
            # if g.inside([polygon_ext], g.Polygon(polygon_int))[0]:
            if g.inside([polygon_int], g.Polygon(polygon_ext))[0]:
                polygons_to_cut.append(polygon_int)
                used_polygons.add(make_hashable(polygon_int))
        if len(polygons_to_cut) >= 2:
            used_polygons.add(
                make_hashable(polygons_to_cut[0])
            )  # don't add an un-cut exterior
            top_level_polygons.append(polygons_to_cut)
    for polygon in polygons:
        if make_hashable(polygon) not in used_polygons:
            top_level_polygons.append([polygon])
    result = [
        reduce(partial(cut_out_hole, layer=layer, datatype=datatype), polygons_to_cut)
        for polygons_to_cut in top_level_polygons
    ]
    return g.PolygonSet(result, layer=layer, datatype=datatype)


@memoize
def character(
    char,
    face=DEFAULT_FACE,
    scale_factor=1,
    points_per_segment=POINTS_PER_SEGMENT,
    layer=0,
    datatype=0,
):
    if char.isalnum():
        char_name = char
    else:
        char_name = "U+{:06x}".format(ord(char))
    cell_name = "~{}~l{}".format(char_name, layer)
    if scale_factor != 1:
        cell_name += "s{}".format(scale_factor)
    cell = Cell(cell_name)
    polygons, metrics = render_character(
        char, face=face, points_per_segment=points_per_segment
    )
    polygons = [p * scale_factor for p in polygons]
    cell.add(cut_out_holes(polygons, layer=layer, datatype=datatype))
    return cell, metrics


def Text(
    text,
    size,
    position=(0, 0),
    horizontal=True,
    angle=0,
    layer=0,
    datatype=0,
    face=DEFAULT_FACE,
    scale_factor=1,
    points_per_segment=POINTS_PER_SEGMENT,
):
    if not horizontal or angle != 0:
        raise NotImplementedError
    position = np.array(position)
    offset = np.zeros(2)
    magnification = size / (face.units_per_EM * scale_factor)
    linespace = face.height * magnification
    face.load_char(" ", freetype.FT_LOAD_NO_SCALE)
    space_advance = face.glyph.metrics.horiAdvance * magnification
    cells = []
    for char in text:
        if char == "\n":
            offset[0] = 0
            offset[1] -= linespace
            continue
        elif char == " ":
            # face.load_char(char)
            offset[0] += space_advance
            continue
        char_cell, metrics = character(
            char,
            face=face,
            scale_factor=scale_factor,
            points_per_segment=points_per_segment,
            layer=layer,
            datatype=datatype,
        )
        cells.append(
            g.CellReference(
                char_cell, position + offset, rotation=0, magnification=magnification
            )
        )
        offset[0] += metrics["horiAdvance"] * magnification
    return cells
