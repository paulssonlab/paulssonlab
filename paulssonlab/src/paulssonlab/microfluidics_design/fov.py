import numpy as np
import pandas as pd
import scipy.optimize as optimize
from itertools import count, cycle, islice, product, accumulate, takewhile
import operator
import networkx as nx
import gdstk
from paulssonlab.microfluidics_design.geometry import Cell
from paulssonlab.microfluidics_design.design import text
from paulssonlab.microfluidics_design.util import get_uuid, memoize
from paulssonlab.util import sign

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
            get_fov_positions(
                fov_grid["origin_y"],
                fov_grid["offsets_y"],
            ),
            num,
        )
    )


def get_fov_positions(origin, offsets, max=None):
    xs = accumulate(
        cycle(offsets),
        operator.add,
        initial=origin,
    )
    if max is not None:
        xs = takewhile(lambda x: x * sign(offsets[0]) <= max * sign(offsets[0]), xs)
    return xs


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


def get_fov_grids_y(
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
            y_prev = None
            offsets_y = []
            permuted_cycle = cycle_[start_idx:] + cycle_[:start_idx]
            margin = np.array([packings_map[y]["total_margin"] for y in permuted_cycle])
            num_active = np.array(
                [packings_map[y]["num_active"] for y in permuted_cycle]
            )
            for idx, y in enumerate(
                islice(cycle(cycle_), start_idx, start_idx + len(cycle_) + 1)
            ):
                if idx == 0:
                    origin_y = y
                    if center_margins:
                        origin_y += packings_map[y]["total_margin"] / 2
                else:
                    offset_y = (
                        packings_map[y_prev]["bottom_interior_active"]
                        - y_prev
                        + packings_map[y]["top_exterior_margin"]
                    )
                    if center_margins:
                        offset_y += (
                            packings_map[y]["total_margin"] / 2
                            - packings_map[y_prev]["total_margin"] / 2
                        )
                    offsets_y.append(offset_y)
                y_prev = y
            fov_grids.append(
                {
                    "origin_y": origin_y,
                    "offsets_y": offsets_y,
                    "margin": margin,
                    "num_active": num_active,
                }
            )
    return fov_grids


def get_fov_grids(fov_dims, metadata, center_margins=True, fov_overlap=None, skip=None):
    if fov_overlap is None:
        fov_overlap = np.array([0, 0])
    else:
        if np.isscalar(fov_overlap):
            fov_overlap = np.array([fov_overlap, fov_overlap])
    grid_metadata = {}
    for fov_name, fov_dim in fov_dims.items():
        grid_metadata[fov_name] = {}
        for region_name, trench_md in metadata["trench_info"].items():
            region_heights = trench_md["region_heights"]
            fov_grids = get_fov_grids_y(
                fov_dim[1],
                trench_md["region_heights"],
                center_margins=center_margins,
                skip=skip,
            )
            # TODO: fov_overlap
            grid_metadata[fov_name][region_name] = []
            for fov_grid in fov_grids:
                offsets_x = [fov_dim[0] - fov_overlap[0]]
                offsets_y = -np.array(fov_grid["offsets_y"])
                ul = trench_md["ul"] + np.array([0, -fov_grid["origin_y"]])
                lr = trench_md["lr"]  # - fov_dim
                xs = get_fov_positions(ul[0], offsets_x, max=lr[0])
                ys = get_fov_positions(ul[1], offsets_y, max=lr[1])
                columns = len(list(xs))
                rows = len(list(ys))
                grid_dim = np.array([columns, rows])
                margin = fov_grid["margin"]
                min_margin = margin.min()
                margin_frac = min_margin / fov_dim[1]
                angle_tol = np.rad2deg(rectangle_rotation_angle(fov_dim, min_margin))
                trench_sets = fov_grid["num_active"]
                grid_metadata[fov_name][region_name].append(
                    {
                        "fov_name": fov_name,
                        "fov_width": fov_dim[0],
                        "fov_height": fov_dim[1],
                        "rows": grid_dim[1],
                        "columns": grid_dim[0],
                        "num_fovs": np.product(grid_dim),
                        "trench_sets": trench_sets,
                        "mean_trench_sets": trench_sets.mean(),
                        "ul": ul,
                        "lr": lr,
                        "offsets_x": offsets_x,
                        "offsets_y": offsets_y,
                        "margin": margin,
                        "min_margin": min_margin,
                        "margin_frac": margin_frac,
                        "angle_tol": angle_tol,  # NOTE: this is in degrees
                        "skip": skip,
                    }
                )
    return grid_metadata


def get_fov_grids_df(
    fov_dims, metadata, center_margins=True, fov_overlap=None, skip=None
):
    grids = get_fov_grids(
        fov_dims,
        metadata,
        center_margins=center_margins,
        fov_overlap=fov_overlap,
        skip=skip,
    )
    return pd.concat(
        {
            fov_name: pd.concat(
                {
                    region: pd.DataFrame(
                        region_grids,
                    ).rename_axis("grid_variant")
                    for region, region_grids in fov_grids.items()
                },
                names=["region"],
            )
            for fov_name, fov_grids in grids.items()
        },
        names=["fov_name"],
    )


@memoize
def draw_fov_cell(
    name, width, height, center_margins=True, angle=None, layer=FOV_LAYER
):
    fov_rect = gdstk.rectangle((0, 0), (width, -height), layer=layer)
    if center_margins:
        rotation_center = (width / 2, -height / 2)
    else:
        rotation_center = (0, 0)
    if angle:
        fov_rect = fov_rect.rotate(angle, rotation_center)
    fov_cell = Cell(f"FOV-{name}")
    fov_cell.add(fov_rect)
    return fov_cell


def draw_fov_grids(
    chip_cell,
    chip_metadata,
    fov_grids_df,
    center_margins=True,
    rotate=False,
    label=True,
    label_font_size=120,
    label_margin=200,
    layer=FOV_LAYER,
):
    fov_layer = layer
    for _, region_fov_grids_df in fov_grids_df.groupby("region"):
        ul = None
        for idx, (key, fov_grid) in enumerate(region_fov_grids_df.iterrows()):
            if ul is None:
                ul = fov_grid["ul"]
            if label:
                label_text = " ".join(
                    f"{k}:{v}" for k, v in zip(fov_grids_df.index.names, key)
                )
                label_anchor = ul + np.array([-label_margin, -idx * label_font_size])
                chip_cell.add(
                    *text(
                        label_text,
                        label_font_size,
                        label_anchor,
                        horizontal_alignment="right",
                        layer=fov_layer,
                    )
                )
            xs = list(
                get_fov_positions(
                    fov_grid["ul"][0], fov_grid["offsets_x"], max=fov_grid["lr"][0]
                )
            )
            ys = list(
                get_fov_positions(
                    fov_grid["ul"][1], fov_grid["offsets_y"], max=fov_grid["lr"][1]
                )
            )
            if rotate:
                angle = np.deg2rad(fov_grid["angle_tol"])
            else:
                angle = None
            fov_cell = draw_fov_cell(
                fov_grid["fov_name"],
                fov_grid["fov_width"],
                fov_grid["fov_height"],
                center_margins=center_margins,
                angle=angle,
                layer=fov_layer,
            )
            for y in ys:
                for x in xs:
                    chip_cell.add(
                        gdstk.Reference(
                            fov_cell,
                            (
                                x,
                                y,
                            ),
                        )
                    )
            fov_layer += 1


def draw_fov_grid_old(
    chip_cell,
    chip_metadata,
    fov_dims,
    grid_metadata,
    rotate=False,
    layer=FOV_LAYER,
):
    for fov_name, fov_layer in zip(fov_dims.keys(), count(layer)):
        fov_dim = fov_dims[fov_name]
        grid_metadata_fov = grid_metadata[fov_name]
        fov_rect = gdstk.rectangle((0, 0), (fov_dim[0], -fov_dim[1]), layer=fov_layer)
        # TODO
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
