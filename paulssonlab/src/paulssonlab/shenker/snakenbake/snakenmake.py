import numpy as np
import gdspy as g
import click
from functools import partial
import numbers

DEFAULT_DIMS = (23e3, 13e3)

Cell = partial(g.Cell, exclude_from_current=True)

# TODO: encode parameters in cell names
# TODO: break out func for each cell, use memoize decorator to reuse cells from multiple chips


def outline(dims, thickness=0.15e3):
    outline_inner = g.Rectangle(-dims / 2, dims / 2)
    outline_outer = g.Rectangle(-(dims + thickness) / 2, (dims + thickness) / 2)
    outline = g.fast_boolean(outline_outer, outline_inner, "not")
    return outline


# def chip(design_func, dims=DEFAULT_DIMS, name=None):
#     cell = Cell('Chip_{}'.format(name) if name else 'Chip')

# TODO: add snake gap using multiple snake calls, put gap between them on chip


def make_odd(number):
    if number % 2 == 0:
        return number - 1
    else:
        return number


def snake(
    dims=DEFAULT_DIMS,
    split=1,
    margin=1e3,
    port_margin=None,
    trench_width=1.5,
    trench_length=35,
    lane_gap=20,
    trench_spacing=2,
    feeding_channel_width=90,
    port_radius=200,
):
    if port_margin is None:
        port_margin = margin / 2
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
    port_cell = Cell("Feeding Channel Port")
    trench_cell = Cell("Trench")
    lane_cell = Cell("Lane")
    snake_cell = Cell("Snake")
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
    return (snake_cell, lane_cell, bend_cell, trench_cell, port_cell), metadata


@click.group()
def cli():
    pass


if __name__ == "__main__":
    cli()
