import numpy as np
from gdspy import Cell, Round, Polygon, boolean
from functools import partial

MAX_POINTS = 4094  # 8191 # same as LayoutEditor
ROUND_POINTS = 100

Cell = partial(Cell, exclude_from_current=True)
Round = partial(Round, number_of_points=ROUND_POINTS, max_points=MAX_POINTS)
boolean = partial(boolean, max_points=MAX_POINTS)


def flatten_or_merge(cell, flatten=False, merge=False, layer=None):
    if flatten or merge:
        cell.flatten()
    if merge:
        cell.polygons = [boolean(cell.polygons, None, "or", layer=layer)]
    return cell


def cross(width, length, **kwargs):
    half_width = width / 2
    points = [
        (half_width, half_width),
        (length, half_width),
        (length, -half_width),
        (half_width, -half_width),
        (half_width, -length),
        (-half_width, -length),
        (-half_width, -half_width),
        (-length, -half_width),
        (-length, half_width),
        (-half_width, half_width),
        (-half_width, length),
        (half_width, length),
    ]
    return Polygon(points, **kwargs)


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
    all_points = np.concatenate(polygons).transpose()
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
            ref.rotation += 180
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
    elif alignment == "centered":
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
    polygons = sum([ref.get_polygons() for ref in refs], [])
    bbox = get_bounding_box(polygons)
    if alignment == "right":
        offset = -np.array([bbox[1, 0], 0])
    elif alignment == "centered":
        offset = -np.array([(bbox[0, 0] + bbox[1, 0]) / 2, 0])
    elif alignment == "left":
        offset = -np.array([bbox[0, 0], 0])
    else:
        raise NotImplementedError
    # position = np.array(position) + offset
    # polygons[:] = [position + factor*poly for poly in polygons]
    # return polygons
    for ref in refs:
        ref.origin += position + offset
    return refs
