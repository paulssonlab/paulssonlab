import numpy as np
import gdspy as g
import click
from functools import partial
from itertools import product
import numbers
from util import make_odd, memoize

TRENCH_LAYER = 2
FEEDING_CHANNEL_LAYER = 6
DEFAULT_DIMS = np.array([23e3, 13e3])

# TODO: freeze numpy arrays, hash dicts: hash(frozenset(self.items()))
Cell = partial(g.Cell, exclude_from_current=True)

# TODO: encode parameters in cell names
# TODO: break out func for each cell, use memoize decorator to reuse cells from multiple chips
# SAMPLER_WAFER FUNCTION TAKES KWARGS, when arrays are given, make variations and name approrpiately

# TODO
# alignment marks
# trenches
# label profilometry marks
# trench numbers (every 100 with tick mark)??
# layers
# union


@memoize
def snake(
    dims=DEFAULT_DIMS,
    split=1,
    margin=1e3,
    port_margin=None,
    trench_width=1.5,
    trench_length=35,
    trench_margin=None,
    lane_gap=20,
    trench_spacing=2,
    feeding_channel_width=90,
    port_radius=200,
    tick_length=5,
    tick_margin=5,
    tick_period=50,
    tick_text_size=None,
    tick_labels=True,
    feeding_channel_layer=FEEDING_CHANNEL_LAYER,
    trench_layer=TRENCH_LAYER,
):
    if port_margin is None:
        port_margin = margin / 2
    if trench_margin is None:
        trench_margin = 10  # min(trench_length, feeding_channel_width / 2)
    if tick_text_size is None:
        tick_text_size = tick_length
    effective_trench_length = trench_length + lane_gap / 2
    inner_snake_turn_radius = effective_trench_length
    outer_snake_turn_radius = feeding_channel_width + inner_snake_turn_radius
    lane_fc_dims = np.array(
        [dims[0] - 2 * margin - outer_snake_turn_radius, feeding_channel_width]
    )
    lane_height = feeding_channel_width + 2 * effective_trench_length
    max_lanes = int((dims[1] - 2 * margin) // lane_height)
    if split is None:
        split = 1
    if isinstance(split, numbers.Integral):
        lanes_per_snake = make_odd(max_lanes // split)
        remainder_lanes_per_snake = make_odd(max_lanes - (split - 1) * lanes_per_snake)
        split = (lanes_per_snake,) * (split - 1) + (remainder_lanes_per_snake,)
    else:
        if isinstance(split[0], numbers.Integral):
            if sum(split) > max_lanes:
                raise Exception(
                    "total lanes desired {} is greater than maximum number of lanes {}".format(
                        sum(split), max_lanes
                    )
                )
            if np.any(np.array(split) % 2 == 0):
                raise Exception("number of lanes per snake must be odd")
        else:
            raise Exception(
                "snake splitting spec must be a integer or sequence of integers"
            )
    num_lanes = sum(split)
    trench_xs = np.arange(
        -(lane_fc_dims[0] - trench_width) / 2,
        (lane_fc_dims[0] + trench_width) / 2,
        trench_width + trench_spacing,
    )
    trenches_per_set = len(trench_xs)
    num_trenches = trenches_per_set * 2 * num_lanes
    metadata = {
        k: v
        for k, v in locals().items()
        if k in ("num_lanes", "trenches_per_set", "num_trenches")
    }
    last_lane_y = ((num_lanes - 1) * lane_height) / 2
    snake_cell = Cell("Snake")
    lane_cell = Cell("Lane")
    lane_fc = g.Rectangle(-lane_fc_dims / 2, lane_fc_dims / 2)
    lane_cell.add(lane_fc)
    bend_cell = Cell("Feeding Channel Bend")
    bend = g.Round(
        (0, 0),
        outer_snake_turn_radius,
        inner_radius=inner_snake_turn_radius,
        number_of_points=500,
    )
    bend = g.slice(bend, 0, 0)
    bend_cell.add(bend[1])
    port_cell = Cell("Feeding Channel Port")
    port_x = dims[0] / 2 - lane_fc_dims[0] / 2 - port_radius - port_margin
    port = g.Round((port_x, 0), port_radius, number_of_points=500)
    port_fc = g.Rectangle((0, -lane_fc_dims[1] / 2), (port_x, lane_fc_dims[1] / 2))
    port_cell.add(port)
    port_cell.add(port_fc)
    lane_ys = np.linspace(last_lane_y, -last_lane_y, num_lanes)
    split_cum = np.concatenate(((0,), np.cumsum(split)))
    right_port_lanes = split_cum[1:] - 1
    left_port_lanes = split_cum[:-1]
    right_bend_lanes = np.concatenate(
        [
            np.arange(start, stop - 1, 2)
            for start, stop in zip(split_cum[:-1], split_cum[1:])
        ]
    )
    left_bend_lanes = right_bend_lanes + 1
    for y in lane_ys:
        snake_cell.add(g.CellReference(lane_cell, (0, y)))
    for lane in right_bend_lanes:
        snake_cell.add(
            g.CellReference(
                bend_cell, (lane_fc_dims[0] / 2, lane_ys[lane] - lane_height / 2)
            )
        )
    for lane in left_bend_lanes:
        snake_cell.add(
            g.CellReference(
                bend_cell,
                (-lane_fc_dims[0] / 2, lane_ys[lane] - lane_height / 2),
                rotation=180,
            )
        )
    for lane in right_port_lanes:
        snake_cell.add(g.CellReference(port_cell, (lane_fc_dims[0] / 2, lane_ys[lane])))
    for lane in left_port_lanes:
        snake_cell.add(
            g.CellReference(
                port_cell, (-lane_fc_dims[0] / 2, lane_ys[lane]), rotation=180
            )
        )
    trench_cell = Cell("Trench")
    trench_cell.add(
        g.Rectangle(
            (-trench_width / 2, -trench_margin),
            (trench_width / 2, trench_length),
            layer=trench_layer,
        )
    )
    tick_cell = Cell("Tick")
    tick_cell.add(
        g.Rectangle(
            (-trench_width / 2, trench_length + tick_margin),
            (trench_width / 2, trench_length + tick_margin + tick_length),
            layer=trench_layer,
        )
    )
    tick_xs = trench_xs[::tick_period]
    num_ticks = len(tick_xs)
    for lane_idx, y in enumerate(lane_ys):
        snake_cell.add(
            g.CellArray(
                trench_cell,
                trenches_per_set,
                1,
                (trench_width + trench_spacing, 0),
                (trench_xs[0], y + feeding_channel_width / 2),
            )
        )
        snake_cell.add(
            g.CellArray(
                trench_cell,
                trenches_per_set,
                1,
                (trench_width + trench_spacing, 0),
                (trench_xs[0], y - feeding_channel_width / 2),
                x_reflection=True,
            )
        )
        snake_cell.add(
            g.CellArray(
                tick_cell,
                num_ticks,
                1,
                (tick_period * (trench_width + trench_spacing), 0),
                (trench_xs[0], y + feeding_channel_width / 2),
            )
        )
        if tick_labels:
            for tick_idx, x in enumerate(tick_xs):
                tick_idx = (
                    lane_idx * 2 * trenches_per_set + 2 * tick_idx * tick_period + 1
                )
                snake_cell.add(
                    g.Text(
                        str(tick_idx),
                        tick_text_size,
                        (
                            x + 2 * trench_width,
                            y + feeding_channel_width / 2 + trench_length + tick_margin,
                        ),
                        layer=trench_layer,
                    )
                )
        # snake_cell.add(g.CellArray(tick_cell, int(trenches_per_set // tick_period), 1, (tick_period*(trench_width + trench_spacing), 0),
        #                            (trench_xs[0], y - feeding_channel_width/2), x_reflection=True))
    return snake_cell, metadata


# TODO: implement using bbox-aware auto-gridding helper
@memoize
def profilometry_marks(
    dims=np.array([150, 75]),
    columns=3,
    rows=3,
    layers=(FEEDING_CHANNEL_LAYER, TRENCH_LAYER),
):
    dims = np.array(dims)
    grid_size = np.array([rows, columns])
    grid_dims = dims * grid_size * 2
    profilometry_cell = Cell("Profilometry Marks")
    mark_cells = {
        layer: Cell("Profilometry Mark Layer {}".format(layer)) for layer in layers
    }
    for i, (layer, cell) in enumerate(mark_cells.items()):
        cell.add(g.Rectangle(dims, -dims / 2, layer=layer))
        origin = -grid_dims / 2 + (i - (len(mark_cells) - 1) / 2) * np.array(
            [grid_dims[0], 0]
        )
        profilometry_cell.add(g.CellArray(cell, columns, rows, dims * 2, origin))
        profilometry_cell.add(
            g.Text(
                str(layer), dims[1], origin - np.array([0, 2 * dims[1]]), layer=layer
            )
        )
    return profilometry_cell


@memoize
def alignment_cross(size=1e3, width=6, layer=TRENCH_LAYER):
    alignment_cell = Cell("Alignment Cross")
    horizontal = g.Rectangle((-size, -width / 2), (size, width / 2))
    vertical = g.Rectangle((-width / 2, -size), (width / 2, size))
    cross = g.fast_boolean(horizontal, vertical, "or", layer=layer)
    alignment_cell.add(cross)
    alignment_cell.add(g.Rectangle((-3 * size, -size), (-2 * size, size), layer=layer))
    alignment_cell.add(g.Rectangle((2 * size, -size), (3 * size, size), layer=layer))
    return alignment_cell


def wafer(
    chips,
    name,
    diameter=76.2e3,
    side=87.15e3,
    chip_area_angle=np.pi / 4,
    feeding_channel_layer=FEEDING_CHANNEL_LAYER,
    trench_layer=TRENCH_LAYER,
):
    if len(chips) > 6:
        raise Exception("cannot lay out more than six chips on a wafer")
    main_cell = Cell("main")
    corner = np.array((side, side)) / 2
    square = g.Rectangle(-corner, corner)
    circle = g.Round((0, 0), diameter / 2, number_of_points=500)
    wafer_outline = g.fast_boolean(square, circle, "not")
    # TODO: put text top horizontal or right vertical depending on chip_angle_area <> np.pi/4
    # TODO: move alignment marks accordingly
    chip_area_corner = (
        diameter / 2 * np.array([np.cos(chip_area_angle), np.sin(chip_area_angle)])
    )
    main_cell.add(g.Rectangle(-chip_area_corner, chip_area_corner))
    main_cell.add(wafer_outline)
    profilometry_cell = profilometry_marks(layers=(feeding_channel_layer, trench_layer))
    profilometry_spacing = np.array([0, chip_area_corner[1] * 2 / 3])
    main_cell.add(
        g.CellArray(
            profilometry_cell, 1, 2, profilometry_spacing, -profilometry_spacing / 2
        )
    )
    alignment_cell = alignment_cross(layer=trench_layer)
    alignment_spacing = np.array([0, chip_area_corner[1] * 7 / 3])
    main_cell.add(
        g.CellArray(alignment_cell, 1, 2, alignment_spacing, -alignment_spacing / 2)
    )
    text_position = (3e4, 1e4)
    main_cell.add(
        g.Text(
            name,
            2000,
            position=text_position,
            angle=-np.pi / 2,
            layer=feeding_channel_layer,
        )
    )
    horizontal_chip_spacing = chip_area_corner[0]
    vertical_chip_spacing = chip_area_corner[1]
    for (horizontal, vertical), chip in zip(product((-1, 1), (-1, 0, 1)), chips):
        x = horizontal_chip_spacing * horizontal / 2
        y = vertical_chip_spacing * vertical * 2 / 3
        main_cell.add(g.CellReference(chip, (x, y)))
    return main_cell


@memoize
def chip(
    name,
    design_func=snake,
    feeding_channel_layer=FEEDING_CHANNEL_LAYER,
    trench_layer=TRENCH_LAYER,
    **kwargs,
):
    if "dims" not in kwargs:
        kwargs["dims"] = DEFAULT_DIMS
    kwargs["feeding_channel_layer"] = feeding_channel_layer
    kwargs["trench_layer"] = trench_layer
    dims = kwargs["dims"]
    design_cell, _ = design_func(**kwargs)
    chip_cell = Cell("Chip-{}".format(name))
    chip_cell.add(outline(dims, layer=feeding_channel_layer))
    chip_cell.add(g.CellReference(design_cell, (0, 0)))
    # text_position = (-1e4,1.2e3)
    text_size = 600
    text_position = (-dims[0] / 2 + 1.5 * text_size, dims[1] / 2 - 1.5 * text_size)
    chip_cell.add(
        g.Text(name, text_size, position=text_position, layer=feeding_channel_layer)
    )
    return chip_cell


@memoize
def outline(dims, thickness=0.15e3, layer=FEEDING_CHANNEL_LAYER):
    outline_inner = g.Rectangle(-dims / 2, dims / 2)
    outline_outer = g.Rectangle(-(dims + thickness) / 2, (dims + thickness) / 2)
    outline = g.fast_boolean(outline_outer, outline_inner, "not", layer=layer)
    return outline


# def chip(design_func, dims=DEFAULT_DIMS, name=None):
#     cell = Cell('Chip_{}'.format(name) if name else 'Chip')

# TODO: add snake gap using multiple snake calls, put gap between them on chip


@click.group()
def cli():
    pass


if __name__ == "__main__":
    cli()
