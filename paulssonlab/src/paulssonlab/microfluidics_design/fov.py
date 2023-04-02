import numpy as np
import scipy.optimize as optimize
from itertools import count, cycle, islice, product, accumulate
import operator
import networkx as nx
import gdstk
from paulssonlab.microfluidics_design.geometry import Cell
from paulssonlab.microfluidics_design.util import get_uuid

FOV_LAYER = 10


def fov_plot(ys, fov_height, region_heights, bar_height=2, bar_spacing=1):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    if np.isscalar(ys):
        ys = [ys]

    total_height = sum(region_heights)
    active_height = len(ys) * bar_height + (len(ys) + 1) * bar_spacing
    fig, ax = plt.subplots(figsize=(16, 4))
    y = 0
    idx = 0
    y_max = np.max(ys) + fov_height
    while y < y_max:
        if idx % 2 == 0:
            ax.add_patch(
                Rectangle(
                    (y, 0),
                    region_heights[idx],
                    active_height,
                    facecolor="tab:blue",
                    lw=0,
                    fill=True,
                )
            )
        y += region_heights[idx]
        idx = (idx + 1) % len(region_heights)
    for fov_num, y in enumerate(reversed(ys)):
        ax.add_patch(
            Rectangle(
                (y, (fov_num + 1) * bar_spacing + fov_num * bar_height),
                fov_height,
                bar_height,
                facecolor="tab:orange",
                fill=True,
            )
        )
    ax.set_xlim((0, y_max))
    ax.set_ylim((0, active_height))


def get_fov_plot_ys(fov_grid, num=4):
    return list(
        islice(
            accumulate(
                cycle(fov_grid["offsets_y"]),
                operator.add,
                initial=fov_grid["origin_y"],
            ),
            num,
        )
    )


def fov_plot_offsets(fov_grids, fov_height, region_heights, num=4, **kwargs):
    if not isinstance(fov_grids, list):
        fov_grids = [fov_grids]
    for fov_grid in fov_grids:
        fov_plot(
            get_fov_plot_ys(fov_grid, num=num), fov_height, region_heights, **kwargs
        )


# TODO: remove
def fov_plot_old(y, fov_height, region_heights):
    import matplotlib.pyplot as plt

    total_height = sum(region_heights)
    repeats = int((2 * fov_height) // total_height)
    plt.figure(figsize=(12, 2))
    plt.step(
        np.concatenate(((0, 0), np.cumsum(region_heights * repeats))),
        np.concatenate(((0, 1), np.arange(len(region_heights) * repeats) % 2 == 0)),
        where="pre",
    )
    plt.plot([0, y, y, y + fov_height, y + fov_height], [0, 0, 0.5, 0.5, 0])


def _next_region(x0, region_heights):
    num_regions = len(region_heights)
    total_height = sum(region_heights)
    if total_height <= 0:
        raise ValueError("total region length must be positive")
    x0_wraparound, x0 = divmod(x0, total_height)
    x = 0
    region = 0
    while True:
        x += region_heights[region]
        if x > x0:
            return x + x0_wraparound * total_height, region
        region = (region + 1) % num_regions


def get_fov_packings(fov_height, region_heights):
    total_height = sum(region_heights)
    packings = []
    x = 0
    while x < total_height:
        new_top, top_idx = _next_region(x, region_heights)
        new_bottom, bottom_idx = _next_region(x + fov_height, region_heights)
        top_interior_margin = new_top - x
        bottom_exterior_margin = new_bottom - (x + fov_height)
        top_exterior_active = new_top - region_heights[top_idx]
        top_exterior_margin = x - top_exterior_active
        top_interior_active = new_top
        bottom_interior_active = new_bottom - region_heights[bottom_idx]
        bottom_interior_margin = x + fov_height - bottom_interior_active
        bottom_exterior_active = new_bottom
        total_margin = min(top_exterior_margin, bottom_interior_margin) + min(
            top_interior_margin, bottom_exterior_margin
        )
        # num_active = (
        #     top_idx - bottom_idx + len(region_heights) * int(fov_height // total_height)
        # )
        num_active = int(
            (top_idx - bottom_idx) % len(region_heights) / 2
            + len(region_heights) / 2 * fov_height // total_height
        )
        if top_idx % 2 == 1 and bottom_idx % 2 == 1:
            packings.append(
                {
                    "top": x,
                    "bottom": x + fov_height,
                    "height": fov_height,
                    "top_exterior_margin": top_exterior_margin,
                    "top_interior_margin": top_interior_margin,
                    "bottom_interior_margin": bottom_interior_margin,
                    "bottom_exterior_margin": bottom_exterior_margin,
                    "top_exterior_active": top_exterior_active,
                    "top_interior_active": top_interior_active,
                    "bottom_exterior_active": bottom_exterior_active,
                    "bottom_interior_active": bottom_interior_active,
                    "total_margin": total_margin,
                    "num_active": num_active,
                }
            )
        if top_interior_margin < bottom_exterior_margin:
            x = new_top
        else:
            x = new_bottom - fov_height
    return packings


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


def _get_fov_graph(packings, total_height):
    graph = nx.DiGraph()
    for packing1, packing2 in product(packings, packings):
        if (
            packing1["bottom_interior_active"] % total_height
            == packing2["top_exterior_active"] % total_height
        ):
            graph.add_edge(packing1["top"], packing2["top"])
    return graph


def get_fov_grids(
    fov_height,
    region_heights,
    center_margins=True,
    skip=None,
):
    if skip is None:
        skip = np.array([0, 0])
    else:
        if np.isscalar(skip):
            skip = np.array([skip, 0])
    total_height = sum(region_heights)
    # TODO: skip is redundant? NO WE NEED TO NUDGE FOV ORIGIN UP
    packings = get_fov_packings(fov_height, region_heights)
    packings_map = {packing["top"]: packing for packing in packings}
    graph = _get_fov_graph(packings, total_height)
    fov_grids = []
    for cycle_ in nx.simple_cycles(graph):
        for start_idx in range(len(cycle_)):
            # offsets_y = []
            # if center_margins:
            #     margin_offset = packings_map[cycle_[start_idx]]["total_margin"] / 2
            # else:
            #     margin_offset = 0
            # origin_y = (
            #     packings_map[cycle_[start_idx]]["top_interior_active"]
            #     - margin_offset % total_height
            # )
            # for y in [
            #     *cycle_[(start_idx + 1) % len(cycle_) :],
            #     *cycle_[: (start_idx + 1) % len(cycle_)],
            # ]:
            # y = origin_y
            y_prev = None
            offsets_y = []
            for idx, y in enumerate(
                islice(cycle(cycle_), start_idx, start_idx + len(cycle_) + 1)
            ):
                if idx == 0:
                    origin_y = y + packings_map[y]["total_margin"] / 2
                else:
                    # offset_y = (
                    #     packings_map[y_prev]["bottom_interior_active"]
                    #     - y_prev
                    #     + packings_map[y]["top_exterior_margin"]
                    #     - packings_map[y_prev]["total_margin"] / 2
                    #     + packings_map[y]["total_margin"] / 2
                    # )
                    offset_y = (
                        packings_map[y_prev]["bottom_interior_active"]
                        - y_prev
                        # - y_prev
                        # + packings_map[y]["top_exterior_active"]
                        + packings_map[y]["top_exterior_margin"]
                        # + y
                        - packings_map[y_prev]["total_margin"] / 2
                        + packings_map[y]["total_margin"] / 2
                        # - packings_map[y]["total_margin"] / 2
                    )
                    offsets_y.append(offset_y)
                # if center_margins:
                #     margin_offset = packings_map[cycle_[start_idx]]["total_margin"] / 2
                # else:
                #     margin_offset = 0
                # margin_offset = 0
                # new_y = (
                #     packings_map[y]["top_interior_active"]
                #     - margin_offset % total_height
                # )
                # offset_y = new_y - old_y
                # offsets_y.append(offset_y)
                # old_y = new_y
                y_prev = y
            fov_grids.append(
                {
                    "origin_y": origin_y,
                    "offsets_y": offsets_y,
                }
            )
    return fov_grids


def get_fov_grids_df(fov_dims, metadata, fov_overlap=None, skip=None):
    if fov_overlap is None:
        fov_overlap = np.array([0, 0])
    else:
        if np.isscalar(fov_overlap):
            fov_overlap = np.array([fov_overlap, fov_overlap])
    grid_metadata = {}
    for fov_name, fov_dim in fov_dims.items():
        for region_name, trench_md in metadata["trench_info"].items():
            region_heights = trench_md["region_heights"]
            fov_grids = get_fov_grids(
                fov_dim,
                trench_md["region_heights"],
                trench_md["origin"],
                trench_md["trench_span"],
                fov_overlap=fov_overlap,
                skip=skip,
            )
            for fov_grid in fov_grids:
                # grid_height = int(np.ceil(metadata["num_lanes"] * 2 / trench_sets_per_fov))
                # max_filled_height = np.max(filled_heights)
                # margin = fov_dim[1] - max_filled_height
                # margin_frac = margin / fov_dim[1]
                offsets_x = [fov_dim[0] - fov_overlap[0]]
                angle_tol = np.rad2deg(rectangle_rotation_angle(fov_dim, margin))
                num_trenches_per_fov = 0
                grid_metadata[fov_name][region_name] = {
                    "fov_name": fov_name,
                    "fov_width": fov_dim[0],
                    "fov_height": fov_dim[1],
                    "fov_area": np.product(fov_dim),
                    "rows": grid_dim[1],
                    "columns": grid_dim[0],
                    "num_fovs": np.product(grid_dim),
                    # "trench_sets_per_fov": trench_sets_per_fov,
                    # "offsets": offsets,
                    # "filled_heights": filled_heights,
                    # "max_filled_heights": max_filled_height,
                    # "margin": margin,
                    # "margin_frac": margin_frac,
                    "angle_tol": angle_tol,  # NOTE: this is in degrees
                    "skip_first": skip_first,
                }
    return grid_metadata


def draw_fov_grid(
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
        ###########
        origin_x = input_metadata["fov_origin_x"]
        origin_y = input_metadata["fov_origin_y"]
        margin = fov_dim[1] - filled_height
        if center_margins:
            margin_offset = margin / 2
        else:
            margin_offset = 0
        y = origin_y + margin_offset
        for offset_y in islice(cycle(offsets_y), num_rows):
            x = origin_x
            for offset_x in islice(cycle(offsets_x), num_columns):
                chip_cell.add(
                    gdstk.Reference(
                        fov_cell,
                        (
                            x,
                            y,
                        ),
                    )
                )
                x += offset_x
            y += offset_y
        ###########
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
