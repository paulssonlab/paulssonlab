import numpy as np
from gdspy import Cell, Round, fast_boolean
from functools import partial

MAX_POINTS = 4094  # 8191 # same as LayoutEditor
ROUND_POINTS = 100

Cell = partial(Cell, exclude_from_current=True)
Round = partial(Round, number_of_points=ROUND_POINTS, max_points=MAX_POINTS)
fast_boolean = partial(fast_boolean, max_points=MAX_POINTS)


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
