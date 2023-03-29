import numpy as np
import scipy.optimize as optimize
from itertools import count
import gdstk
from paulssonlab.shenker.snakenbake.geometry import Cell
from paulssonlab.shenker.snakenbake.util import get_uuid

FOV_LAYER = 10


def get_fov_coverage(
    fov,
    feeding_channel_width=None,
    trench_length=None,
    trench_gap=None,
    num_lanes=None,
    lane_with_trenches_length=None,
    trench_sets_per_fov=np.inf,
    skip_first=False,
    **kwargs,
):
    trench_sets_per_fov_target = trench_sets_per_fov
    trench_sets_per_fov = 1
    start_by_skipping = skip_first is True
    y = trench_length
    unit_cell_height = None
    while True:
        delta_y = trench_gap + trench_length
        if y + delta_y < fov[1] and trench_sets_per_fov < trench_sets_per_fov_target:
            if start_by_skipping is not False:
                start_by_skipping = None
                y += delta_y
                offset = y + feeding_channel_width
                trench_sets_per_fov += 1
        else:
            break
        delta_y = feeding_channel_width + trench_length
        if y + delta_y < fov[1] and trench_sets_per_fov < trench_sets_per_fov_target:
            if start_by_skipping is not True:
                start_by_skipping = None
                y += delta_y
                offset = y + trench_gap
                trench_sets_per_fov += 1
        else:
            break
    return trench_sets_per_fov, offset, y


def rectangle_rotation_angle(fov_dim, margin):
    width, height = fov_dim
    return 2 * np.arctan(
        (np.sqrt(2 * height * margin + width**2 - margin**2) - width)
        / (2 * height - margin)
    )


def _rectangle_rotation_angle_root_finder(fov_dim, margin):
    width, height = fov_dim

    def func(theta):
        return width * np.sin(theta) - height * (np.cos(theta) - 1) - margin

    res = optimize.root_scalar(func, bracket=(0, np.pi / 4))
    if not res.converged:
        return None
    else:
        return res.root


def get_grid_metadata(fov_dims, metadata, skip_first=None):
    grid_metadata = {}
    for fov_name, fov_dim in fov_dims.items():
        trench_sets_per_fov, offset, filled_height = get_fov_coverage(
            fov_dim, **metadata, skip_first=skip_first
        )
        if skip_first is None:
            trench_sets_per_fov2, offset2, filled_height2 = get_fov_coverage(
                fov_dim, **metadata, skip_first=skip_first
            )
            # if skip_first=None, pick True/False depending on which gives more trench_sets_per_fov
            if trench_sets_per_fov2 > trench_sets_per_fov:
                trench_sets_per_fov, offset, filled_height = (
                    trench_sets_per_fov2,
                    offset2,
                    filled_height2,
                )
        if trench_sets_per_fov % 2 != 0:
            # TODO: there should be a more elegant way of getting this without a second call
            _, offset2, filled_height2 = get_fov_coverage(
                fov_dim,
                **metadata,
                trench_sets_per_fov=trench_sets_per_fov,
                skip_first=not skip_first,
            )
            offsets = (offset, offset2)
            filled_heights = (filled_height, filled_height2)
        else:
            offsets = (offset,)
            filled_heights = (filled_height,)
        grid_width = int(np.ceil(metadata["lane_with_trenches_length"] / fov_dim[0]))
        grid_height = int(np.ceil(metadata["num_lanes"] * 2 / trench_sets_per_fov))
        max_filled_height = np.max(filled_heights)
        margin = fov_dim[1] - max_filled_height
        margin_frac = margin / fov_dim[1]
        angle_tol = np.rad2deg(rectangle_rotation_angle(fov_dim, margin))
        grid_metadata[fov_name] = {
            "fov_name": fov_name,
            "fov_width": fov_dim[0],
            "fov_height": fov_dim[1],
            "fov_area": np.product(fov_dim),
            "grid_width": grid_width,
            "grid_height": grid_height,
            "num_fovs": grid_width * grid_height,
            "trench_sets_per_fov": trench_sets_per_fov,
            "offsets": offsets,
            "filled_heights": filled_heights,
            "max_filled_heights": max_filled_height,
            "margin": margin,
            "margin_frac": margin_frac,
            "angle_tol": angle_tol,  # NOTE: this is in degrees
            "skip_first": skip_first,
        }
    return grid_metadata


def draw_grid_overlay(
    chip_cell,
    chip_metadata,
    fov_dims,
    grid_metadata,
    center_margins=True,
    rotate=False,
    layer=FOV_LAYER,
):
    for fov_name, fov_layer in zip(fov_dims.keys(), count(layer)):
        fov_dim = fov_dims[fov_name]
        grid_metadata_fov = grid_metadata[fov_name]
        fov_rect = gdstk.rectangle((0, 0), (fov_dim[0], -fov_dim[1]), layer=fov_layer)
        if center_margins:
            rotation_center = (fov_dim[0] / 2, -fov_dim[1] / 2)
        else:
            rotation_center = (0, 0)
        if rotate:
            fov_rect = fov_rect.rotate(
                np.deg2rad(grid_metadata_fov["angle_tol"]), rotation_center
            )
        fov_cell = Cell(f"FOV-{fov_name}")
        fov_cell.add(fov_rect)
        y = 0
        offsets = grid_metadata_fov["offsets"]
        filled_heights = grid_metadata_fov["filled_heights"]
        total_offset = np.sum(offsets)
        for idx, (offset, filled_height) in enumerate(zip(offsets, filled_heights)):
            margin = fov_dim[1] - filled_height
            if center_margins:
                margin_offset = margin / 2
            else:
                margin_offset = 0
            chip_cell.add(
                gdstk.Reference(
                    fov_cell,
                    (
                        chip_metadata["fov_origin_x"],
                        chip_metadata["fov_origin_y"] + margin_offset + y,
                    ),
                    columns=grid_metadata_fov["grid_width"],
                    rows=(grid_metadata_fov["grid_height"] - (idx + 1)) // len(offsets)
                    + 1,
                    spacing=(fov_dim[0], -total_offset),
                )
            )
            y -= offset
