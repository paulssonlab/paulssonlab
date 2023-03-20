import numpy as np
import gdstk
from matplotlib.path import Path
from functools import partial
from paulssonlab.shenker.snakenbake.util import get_uuid, get_polygons

ELLIPSE_TOLERANCE = 0.3


def Cell(name):
    if name != "main":
        name = f"{name}-{get_uuid()}"
    return gdstk.Cell(name)


ellipse = partial(gdstk.ellipse, tolerance=ELLIPSE_TOLERANCE)


def flatten_or_merge(cell, flatten=False, merge=False, layer=None):
    if flatten or merge:
        cell.flatten()
    if merge:
        polygons = gdstk.boolean(cell.polygons, [], "or", layer=layer)
        cell.remove(*cell.polygons)
        cell.add(*polygons)
    return cell


def cross(length, thickness, **kwargs):
    half_thickness = thickness / 2
    points = [
        (half_thickness, half_thickness),
        (length, half_thickness),
        (length, -half_thickness),
        (half_thickness, -half_thickness),
        (half_thickness, -length),
        (-half_thickness, -length),
        (-half_thickness, -half_thickness),
        (-length, -half_thickness),
        (-length, half_thickness),
        (-half_thickness, half_thickness),
        (-half_thickness, length),
        (half_thickness, length),
    ]
    return gdstk.Polygon(points, **kwargs)


def l_shape(height, width, thickness, **kwargs):
    half_thickness = thickness / 2
    points = [
        (-half_thickness, height),
        (half_thickness, height),
        (half_thickness, -(height - thickness)),
        (width - half_thickness, -(height - thickness)),
        (width - half_thickness, -height),
        (-half_thickness, -height),
    ]
    return gdstk.Polygon(points, **kwargs)


def qr_target(outer_thickness, margin, inner_width, layer=None):
    hole_x = margin + inner_width / 2
    outer_x = hole_x + outer_thickness
    half_inner_width = inner_width / 2
    outer = gdstk.rectangle((-outer_x, -outer_x), (outer_x, outer_x), layer=layer)
    hole = gdstk.rectangle((-hole_x, -hole_x), (hole_x, hole_x), layer=layer)
    outer = gdstk.boolean(outer, hole, "not", layer=layer)
    inner = gdstk.rectangle(
        (-half_inner_width, -half_inner_width),
        (half_inner_width, half_inner_width),
        layer=layer,
    )
    return gdstk.boolean(outer, inner, "or", layer=layer)


def polygon_orientation(polygon):
    # SEE: https://en.wikipedia.org/wiki/Curve_orientation
    polygon = np.array(polygon)
    x_min = polygon[:, 0].min()
    y_min = polygon[polygon[:, 0] == x_min, 1].min()
    idx_b = np.where((polygon[:, 0] == x_min) & (polygon[:, 1] == y_min))[0][0]
    idx_a = (idx_b - 1) % len(polygon)
    idx_c = (idx_b + 1) % len(polygon)
    x_a, y_a = polygon[idx_a]
    x_b, y_b = polygon[idx_b]
    x_c, y_c = polygon[idx_c]
    det = (x_b - x_a) * (y_c - y_a) - (x_c - x_a) * (y_b - y_c)
    return det


def get_bounding_box(polygons):
    all_points = np.concatenate([p.points for p in polygons]).transpose()
    bbox = np.array(
        (
            (all_points[0].min(), all_points[1].min()),
            (all_points[0].max(), all_points[1].max()),
        )
    )
    return bbox


def mirror(polygons, axis="y"):
    x = y = False
    if axis == "both":
        x = y = True
    elif axis == "x":
        x = True
    elif axis == "y":
        y = True
    else:
        raise ValueError("bad axis specification")
    factor = np.ones(2)
    if x:
        factor[0] = -1
    if y:
        factor[1] = -1
    polygons[:] = [factor * poly for poly in polygons]
    return polygons


def mirror_refs(refs, axis="x"):
    x = y = False
    if axis == "both":
        x = y = True
    elif axis == "x":
        x = True
    elif axis == "y":
        y = True
    else:
        raise ValueError("bad axis specification")
    factor = np.ones(2)
    if x:
        factor[0] = -1
    if y:
        factor[1] = -1
    for ref in refs:
        if x:
            ref.rotation += np.deg2rad(180)
            ref.x_reflection = not ref.x_reflection
        if y:
            raise NotImplementedError
        ref.origin *= factor
    return refs


# TODO: right/left switched
def align(polygons, position=(0, 0), alignment="left"):
    bbox = get_bounding_box(polygons)
    if alignment == "left":
        offset = np.array([bbox[1, 0], 0])
    elif alignment == "center":
        offset = np.array([(bbox[0, 0] + bbox[1, 0]) / 2, 0])
    elif alignment == "right":
        offset = np.array([bbox[0, 0], 0])
    else:
        raise NotImplementedError
    factor = np.array([-1, 1])
    position = np.array(position) + offset
    polygons[:] = [position + factor * poly for poly in polygons]
    return polygons


def align_refs(refs, position=(0, 0), alignment="left"):
    position = np.array(position)
    polygons = get_polygons(refs)
    bbox = get_bounding_box(polygons)
    if alignment == "right":
        offset = -np.array([bbox[1, 0], 0])
    elif alignment == "center":
        offset = -np.array([(bbox[0, 0] + bbox[1, 0]) / 2, 0])
    elif alignment == "left":
        offset = -np.array([bbox[0, 0], 0])
    else:
        raise NotImplementedError
    for ref in refs:
        ref.origin = np.array(ref.origin) + position + offset
    return refs


# FROM: https://heitzmann.github.io/g/how-tos.html#system-fonts
def from_matplotlib_path(path, layer=0, datatype=0, tolerance=0.1):
    precision = 0.1 * tolerance
    x_max = None
    polys = []
    for points, code in path.iter_segments():
        if code == path.MOVETO:
            c = gdstk.Curve(points, tolerance=tolerance)
        elif code == path.LINETO:
            c.segment(points.reshape(points.size // 2, 2))
        elif code == path.CURVE3:
            c.quadratic(points.reshape(points.size // 2, 2))
        elif code == path.CURVE4:
            c.cubic(points.reshape(points.size // 2, 2))
        elif code == path.CLOSEPOLY:
            pts = c.points()
            if pts.size > 0:
                poly = gdstk.Polygon(pts, layer=layer, datatype=datatype)
                if x_max is not None and pts[:, 0].min() < x_max:
                    i = len(polys) - 1
                    while i >= 0:
                        if polys[i].contain_any(*poly.points):
                            p = polys.pop(i)
                            poly = gdstk.boolean(
                                p,
                                poly,
                                "xor",
                                precision,
                                layer=layer,
                                datatype=datatype,
                            )[0]
                            break
                        elif poly.contain_any(*polys[i].points):
                            p = polys.pop(i)
                            poly = gdstk.boolean(
                                p,
                                poly,
                                "xor",
                                precision,
                                layer=layer,
                                datatype=datatype,
                            )[0]
                        i -= 1
                x_max_new = poly.points[:, 0].max()
                if x_max is None:
                    x_max = x_max_new
                else:
                    x_max = max(x_max, x_max_new)
                polys.append(poly)
    return polys


def from_matplotlib_path_codes(verts, codes, tolerance=0.1, layer=0, datatype=0):
    return from_matplotlib_path(
        Path(verts, codes, closed=False),
        tolerance=tolerance,
        layer=layer,
        datatype=datatype,
    )
