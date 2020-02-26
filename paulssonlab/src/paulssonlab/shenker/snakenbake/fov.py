import numpy as np


def fovs_per_chip(
    fov,
    feeding_channel_width=None,
    trench_length=None,
    trench_gap=None,
    num_lanes=None,
    lane_with_trenches_length=None,
    **kwargs,
):
    hfovs = int(np.ceil(lane_with_trenches_length / fov[0]))
    y = 2 * trench_length + trench_gap
    unit_cell_height = y + feeding_channel_width
    trench_sets_per_fov = 2
    while True:
        delta_y = feeding_channel_width + trench_length
        if y + delta_y < fov[1]:
            y += delta_y
            unit_cell_height = y + trench_gap
            trench_sets_per_fov += 1
        else:
            break
        delta_y = trench_gap + trench_length
        if y + delta_y < fov[1]:
            y += delta_y
            unit_cell_height = y + feeding_channel_width
            trench_sets_per_fov += 1
        else:
            break
    vfovs = int(np.ceil(num_lanes * 2 / trench_sets_per_fov))
    return np.array([hfovs, vfovs]), unit_cell_height, y
