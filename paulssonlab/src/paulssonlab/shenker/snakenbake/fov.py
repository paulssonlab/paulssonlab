import numpy as np
from itertools import count
import gdspy
from geometry import Cell

FOV_LAYER = 10


def get_fov_coverage(
    fov,
    feeding_channel_width=None,
    trench_length=None,
    trench_gap=None,
    num_lanes=None,
    lane_with_trenches_length=None,
    skip_first=False,
    **kwargs,
):
    start_by_skipping = skip_first
    y = trench_length
    unit_cell_height = None
    trench_sets_per_fov = 1
    while True:
        delta_y = trench_gap + trench_length
        if y + delta_y < fov[1]:
            if start_by_skipping is not False:
                start_by_skipping = None
                y += delta_y
                offset = y + feeding_channel_width
                trench_sets_per_fov += 1
        else:
            break
        delta_y = feeding_channel_width + trench_length
        if y + delta_y < fov[1]:
            if start_by_skipping is not True:
                start_by_skipping = None
                y += delta_y
                offset = y + trench_gap
                trench_sets_per_fov += 1
        else:
            break
    return trench_sets_per_fov, offset, y


def get_grid_metadata(fov_dims, metadata, skip_first=False):
    grid_metadata = {}
    for fov_name, fov_dim in fov_dims.items():
        trench_sets_per_fov, offset1, active_height = get_fov_coverage(
            fov_dim, **metadata, skip_first=skip_first
        )
        if trench_sets_per_fov % 2 != 0:
            # TODO: there should be a more elegant way of getting this without a second call
            _, offset2, _ = get_fov_coverage(
                fov_dim, **metadata, skip_first=not skip_first
            )
            offsets = (offset1, offset2)
        else:
            offsets = (offset1,)
        grid_width = int(np.ceil(metadata["lane_with_trenches_length"] / fov_dim[0]))
        grid_height = int(np.ceil(metadata["num_lanes"] * 2 / trench_sets_per_fov))
        margin = fov_dim[1] - active_height
        margin_frac = margin / fov_dim[1]
        angle_tol = np.rad2deg(np.arccos(active_height / fov_dim[1]))
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
            "active_height": active_height,
            "margin": margin,
            "margin_frac": margin_frac,
            "angle_tol": angle_tol,
            "skip_first": skip_first,
        }
    return grid_metadata


def draw_grid_overlay(
    chip_cell, chip_metadata, fov_dims, grid_metadata, layer=FOV_LAYER
):
    for fov_name, fov_layer in zip(fov_dims.keys(), count(layer)):
        fov_dim = fov_dims[fov_name]
        grid_metadata_fov = grid_metadata[fov_name]
        fov_rect = gdspy.Rectangle((0, 0), (fov_dim[0], -fov_dim[1]), layer=fov_layer)
        fov_cell = Cell(f"FOV-{fov_name}")
        fov_cell.add(fov_rect)
        y = 0
        offsets = grid_metadata_fov["offsets"]
        total_offset = np.sum(offsets)
        for idx, offset in enumerate(offsets):
            chip_cell.add(
                gdspy.CellArray(
                    fov_cell,
                    grid_metadata_fov["grid_width"],
                    (grid_metadata_fov["grid_height"] - (idx + 1)) // len(offsets) + 1,
                    (fov_dim[0], -total_offset),
                    (chip_metadata["fov_origin_x"], chip_metadata["fov_origin_y"] + y),
                )
            )
            y -= offset
