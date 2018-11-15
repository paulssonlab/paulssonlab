#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import gdspy as g
from gdspy import CellReference, CellArray, Rectangle
from geometry import (
    MAX_POINTS,
    ROUND_POINTS,
    Cell,
    Round,
    fast_boolean,
    mirror,
    mirror_refs,
    align,
    align_refs,
)
from text import Text as _Text
from util import make_odd, memoize
import click
from functools import partial
from cytoolz import compose
from itertools import product
import numbers
import shortuuid

TRENCH_LAYER = 2
FEEDING_CHANNEL_LAYER = 6
DEFAULT_DIMS = np.array([23e3, 13e3])


def Text(
    text,
    size,
    position=(0, 0),
    alignment="left",
    prerotate_alignment=None,
    angle=0,
    **kwargs,
):
    objs = _Text(text, size, position=(0, 0), **kwargs)
    objs = mirror_refs(objs)
    if prerotate_alignment is not None:
        objs = align_refs(objs, position=(0, 0), alignment=prerotate_alignment)
    if angle == 0:
        pass
    elif angle == -np.pi / 2:
        for ref in objs:
            ref.rotation += 90
            ref.origin[:] = ref.origin[::-1]
    elif angle == np.pi / 2:
        for ref in objs:
            ref.rotation -= 90
            ref.origin[:] = ref.origin[::-1] * np.array([1, -1])
    else:
        raise NotImplementedError
    objs = align_refs(objs, position=position, alignment=alignment)
    return objs


get_uuid = partial(shortuuid.random, length=2)

# TODO: encode parameters in cell names
# TODO: break out func for each cell, use memoize decorator to reuse cells from multiple chips
# SAMPLER_WAFER FUNCTION TAKES KWARGS, when arrays are given, make variations and name approrpiately

# TODO
# autonaming of cells according to argument that changed, autolabeling of sampler_wafer
# make ROUND_POINST, MAX_POINTS modifyable at runtime by replacing partial

# .5mm margin outside of port for punching
# .5mm mixing zone after bends
# use round number for align mark value, brandon uses -32000um
# fix font to reduce x-extent
# make index numbers larger so smallest feature is 1.5um across
# make alternate wafer with 3 only chips, centered: centered
# put 1/bottom on bottom align mark
# mirror all text


@memoize
def snake(
    dims=DEFAULT_DIMS,
    split=1,
    horizontal_margin=2.5e3,
    top_margin=1.5e3,
    bottom_margin=1e3,
    port_margin=None,
    trench_width=1.5,
    trench_length=35,
    trench_fc_overlap=None,
    trench_margin=0.5e3,
    lane_gap=20,
    trench_spacing=2,
    feeding_channel_width=90,
    port_radius=200,
    tick_length=5,
    tick_margin=5,
    tick_period=50,
    tick_text_size=None,
    tick_labels=True,
    draw_trenches=True,
    flatten_feeding_channel=False,
    merge_feeding_channel=True,
    feeding_channel_layer=FEEDING_CHANNEL_LAYER,
    trench_layer=TRENCH_LAYER,
    label=None,
):
    if label is None:
        label = get_uuid()
    if port_margin is None:
        port_margin = horizontal_margin / 2
    if trench_fc_overlap is None:
        trench_fc_overlap = min(trench_length, feeding_channel_width / 3)
    if tick_text_size is None:
        tick_text_size = tick_length * 2
    effective_trench_length = trench_length + lane_gap / 2
    inner_snake_turn_radius = effective_trench_length
    outer_snake_turn_radius = feeding_channel_width + inner_snake_turn_radius
    lane_fc_dims = np.array(
        [
            dims[0] - 2 * horizontal_margin - outer_snake_turn_radius,
            feeding_channel_width,
        ]
    )
    lane_height = feeding_channel_width + 2 * effective_trench_length
    max_lanes = int((dims[1] - top_margin - bottom_margin) // lane_height)
    if split is None:
        split = 1
    if isinstance(split, numbers.Integral):
        good_split = False
        for lanes_per_snake in (int(np.around(max_lanes / split)), max_lanes // split):
            lanes_per_snake = make_odd(lanes_per_snake)
            remainder_lanes_per_snake = make_odd(
                max_lanes - (split - 1) * lanes_per_snake
            )
            new_split = (lanes_per_snake,) * (split - 1) + (remainder_lanes_per_snake,)
            if all(s > 1 for s in new_split):
                good_split = True
                split = new_split
                break
        if not good_split:
            raise ValueError("bad split: {split}".format(",".join(split)))
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
        -(lane_fc_dims[0] - trench_margin - trench_width) / 2,
        (lane_fc_dims[0] - trench_margin - trench_width) / 2,
        trench_width + trench_spacing,
    )
    trenches_per_set = len(trench_xs)
    num_trenches = trenches_per_set * 2 * num_lanes
    metadata = {
        k: v
        for k, v in locals().items()
        if k
        in (
            "num_lanes",
            "trenches_per_set",
            "num_trenches",
            "split",
            "feeding_channel_width",
        )
    }
    metadata["lane_length"] = lane_fc_dims[0]
    last_lane_y = ((num_lanes - 1) * lane_height) / 2
    snake_cell = Cell("Snake-{}".format(label))
    lane_cell = Cell("Lane-{}".format(label))
    lane_fc = Rectangle(
        -lane_fc_dims / 2, lane_fc_dims / 2, layer=feeding_channel_layer
    )
    lane_cell.add(lane_fc)
    bend_cell = Cell("Feeding Channel Bend-{}".format(label))
    bend = Round(
        (0, 0),
        outer_snake_turn_radius,
        inner_radius=inner_snake_turn_radius,
        layer=feeding_channel_layer,
    )
    bend = g.slice(bend, 0, 0, layer=feeding_channel_layer)
    bend_cell.add(bend[1])
    port_cell = Cell("Feeding Channel Port-{}".format(label))
    port_x = dims[0] / 2 - lane_fc_dims[0] / 2 - port_radius - port_margin
    port = Round((port_x, 0), port_radius, layer=feeding_channel_layer)
    port_fc = Rectangle(
        (0, -lane_fc_dims[1] / 2),
        (port_x, lane_fc_dims[1] / 2),
        layer=feeding_channel_layer,
    )
    port_cell.add(port)
    port_cell.add(port_fc)
    lane_ys = (
        np.linspace(last_lane_y, -last_lane_y, num_lanes)
        + (bottom_margin - top_margin) / 2
    )
    split_cum = np.concatenate(((0,), np.cumsum(split)))
    left_port_lanes = split_cum[1:] - 1
    right_port_lanes = split_cum[:-1]
    left_bend_lanes = np.concatenate(
        [
            np.arange(start, stop - 1, 2)
            for start, stop in zip(split_cum[:-1], split_cum[1:])
        ]
    )
    right_bend_lanes = left_bend_lanes + 1
    snake_fc_cell = Cell("Snake Feeding Channel-{}".format(label))
    for y in lane_ys:
        snake_fc_cell.add(CellReference(lane_cell, (0, y)))
    for lane in right_bend_lanes:
        snake_fc_cell.add(
            CellReference(
                bend_cell, (lane_fc_dims[0] / 2, lane_ys[lane] - lane_height / 2)
            )
        )
    for lane in left_bend_lanes:
        snake_fc_cell.add(
            CellReference(
                bend_cell,
                (-lane_fc_dims[0] / 2, lane_ys[lane] - lane_height / 2),
                rotation=180,
            )
        )
    for lane in right_port_lanes:
        snake_fc_cell.add(
            CellReference(port_cell, (lane_fc_dims[0] / 2, lane_ys[lane]))
        )
    for lane in left_port_lanes:
        snake_fc_cell.add(
            CellReference(
                port_cell, (-lane_fc_dims[0] / 2, lane_ys[lane]), rotation=180
            )
        )
    if flatten_feeding_channel or merge_feeding_channel:
        snake_fc_cell.flatten()
    if merge_feeding_channel:
        assert len(snake_fc_cell.elements) == 1
        snake_fc_cell.elements[0] = fast_boolean(
            snake_fc_cell.elements[0], None, "or", layer=feeding_channel_layer
        )
        # print(len(snake_fc_cell.elements[0].polygons))
    snake_cell.add(CellReference(snake_fc_cell, (0, 0)))
    trench_cell = Cell("Trench-{}".format(label))
    trench_cell.add(
        Rectangle(
            (-trench_width / 2, -trench_fc_overlap),
            (trench_width / 2, trench_length),
            layer=trench_layer,
        )
    )
    tick_cell = Cell("Tick-{}".format(label))
    tick_cell.add(
        Rectangle(
            (-trench_width / 2, trench_length + tick_margin),
            (trench_width / 2, trench_length + tick_margin + tick_length),
            layer=trench_layer,
        )
    )
    tick_xs = trench_xs[::tick_period]
    num_ticks = len(tick_xs)
    if draw_trenches:
        for lane_idx, y in enumerate(lane_ys):
            snake_cell.add(
                CellArray(
                    trench_cell,
                    trenches_per_set,
                    1,
                    (trench_width + trench_spacing, 0),
                    (trench_xs[0], y + feeding_channel_width / 2),
                )
            )
            snake_cell.add(
                CellArray(
                    trench_cell,
                    trenches_per_set,
                    1,
                    (trench_width + trench_spacing, 0),
                    (trench_xs[0], y - feeding_channel_width / 2),
                    x_reflection=True,
                )
            )
            snake_cell.add(
                CellArray(
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
                        Text(
                            str(tick_idx),
                            tick_text_size,
                            (
                                x + 2 * trench_width,
                                y
                                + feeding_channel_width / 2
                                + trench_length
                                + tick_margin,
                            ),
                            layer=trench_layer,
                        )
                    )
            # snake_cell.add(CellArray(tick_cell, int(trenches_per_set // tick_period), 1, (tick_period*(trench_width + trench_spacing), 0),
            #                            (trench_xs[0], y - feeding_channel_width/2), x_reflection=True))
    return snake_cell, metadata


# TODO: implement using bbox-aware auto-gridding helper
@memoize
def profilometry_marks(
    dims=np.array([150, 75]),
    columns=3,
    rows=3,
    layers=(FEEDING_CHANNEL_LAYER, TRENCH_LAYER),
    text=True,
):
    dims = np.array(dims)
    grid_size = np.array([rows, columns])
    grid_dims = dims * grid_size * 2
    profilometry_cell = Cell("Profilometry Marks")
    mark_cells = {
        layer: Cell("Profilometry Mark Layer {}".format(layer)) for layer in layers
    }
    for i, (layer, cell) in enumerate(mark_cells.items()):
        cell.add(Rectangle(dims, -dims / 2, layer=layer))
        origin = -grid_dims / 2 + (i - (len(mark_cells) - 1) / 2) * np.array(
            [grid_dims[0], 0]
        )
        profilometry_cell.add(CellArray(cell, columns, rows, dims * 2, origin))
        if text:
            profilometry_cell.add(
                Text(
                    str(layer),
                    dims[1],
                    origin - np.array([0, 2 * dims[1]]),
                    layer=layer,
                )
            )
    return profilometry_cell


@memoize
def alignment_cross(size=1e3, width=6, layer=TRENCH_LAYER):
    alignment_cell = Cell("Alignment Cross")
    horizontal = Rectangle((-size, -width / 2), (size, width / 2))
    vertical = Rectangle((-width / 2, -size), (width / 2, size))
    cross = fast_boolean(horizontal, vertical, "or", layer=layer)
    alignment_cell.add(cross)
    alignment_cell.add(Rectangle((-3 * size, -size), (-2 * size, size), layer=layer))
    alignment_cell.add(Rectangle((2 * size, -size), (3 * size, size), layer=layer))
    return alignment_cell


def wafer(
    chips,
    name,
    diameter=76.2e3,
    side=87.15e3,
    chip_area_angle=np.pi / 4,
    chip_area_margin=4e3,
    alignment_mark_position=32e3,
    alignment_text_size=1000,
    label_text_size=2000,
    text=True,
    feeding_channel_layer=FEEDING_CHANNEL_LAYER,
    trench_layer=TRENCH_LAYER,
):
    if len(chips) > 6:
        raise Exception("cannot lay out more than six chips on a wafer")
    main_cell = Cell("main")
    corner = np.array((side, side)) / 2
    square = Rectangle(-corner, corner)
    circle = Round((0, 0), diameter / 2)
    wafer_outline = fast_boolean(square, circle, "not")
    # TODO: put text top horizontal or right vertical depending on chip_angle_area <> np.pi/4
    # TODO: move alignment marks accordingly
    chip_area_corner = (diameter / 2 - chip_area_margin) * np.array(
        [np.cos(chip_area_angle), np.sin(chip_area_angle)]
    )
    if alignment_mark_position is None:
        alignment_mark_position = chip_area_corner[1] * 7 / 6
    main_cell.add(Rectangle(-chip_area_corner, chip_area_corner))
    main_cell.add(wafer_outline)
    profilometry_cell = profilometry_marks(
        layers=(feeding_channel_layer, trench_layer), text=text
    )
    profilometry_spacing = np.array([0, chip_area_corner[1] * 2 / 3])
    main_cell.add(
        CellArray(
            profilometry_cell, 1, 2, profilometry_spacing, -profilometry_spacing / 2
        )
    )
    alignment_cell = alignment_cross(layer=trench_layer)
    alignment_spacing = np.array([0, alignment_mark_position * 2])
    if text:
        main_cell.add(
            Text(
                "bottom",
                alignment_text_size,
                position=(0, 2 * alignment_text_size - alignment_spacing[1] / 2),
                alignment="centered",
                layer=trench_layer,
            )
        )
    main_cell.add(
        CellArray(alignment_cell, 1, 2, alignment_spacing, -alignment_spacing / 2)
    )
    # main_cell.add(Text(name, label_text_size, position=text_position, alignment='centered', angle=-np.pi/2, layer=feeding_channel_layer))
    text_position = (chip_area_corner[0] + label_text_size, 0)
    if text:
        main_cell.add(
            Text(
                name,
                label_text_size,
                position=text_position,
                angle=np.pi / 2,
                alignment="left",
                prerotate_alignment="centered",
                layer=feeding_channel_layer,
            )
        )
    horizontal_chip_spacing = chip_area_corner[0]
    vertical_chip_spacing = chip_area_corner[1]
    if len(chips) <= 3:
        horizontal_spacings = (0,)
    else:
        horizontal_spacings = (-1, 1)
    for (horizontal, vertical), chip in zip(
        product(horizontal_spacings, (-1, 0, 1)), chips
    ):
        x = horizontal_chip_spacing * horizontal / 2
        y = vertical_chip_spacing * vertical * 2 / 3
        main_cell.add(CellReference(chip, (x, y)))
    return main_cell


def chip(
    name,
    design_func=snake,
    label_text_size=600,
    text=True,
    feeding_channel_layer=FEEDING_CHANNEL_LAYER,
    trench_layer=TRENCH_LAYER,
    metadata=None,
    **kwargs,
):
    if "dims" not in kwargs:
        kwargs["dims"] = DEFAULT_DIMS
    kwargs["feeding_channel_layer"] = feeding_channel_layer
    kwargs["trench_layer"] = trench_layer
    dims = kwargs["dims"]
    design_cell, md = design_func(**{"label": name, **kwargs})
    if metadata is not None:
        metadata[name] = md
    chip_cell = Cell("Chip-{}".format(name))
    chip_cell.add(outline(dims, layer=feeding_channel_layer))
    chip_cell.add(CellReference(design_cell, (0, 0)))
    text_position = (0, dims[1] / 2 - 1.5 * label_text_size)
    if text:
        chip_cell.add(
            Text(
                name,
                label_text_size,
                position=text_position,
                alignment="centered",
                layer=feeding_channel_layer,
            )
        )
    return chip_cell


@memoize
def outline(dims, thickness=0.15e3, layer=FEEDING_CHANNEL_LAYER):
    outline_inner = Rectangle(-dims / 2, dims / 2)
    outline_outer = Rectangle(-(dims + thickness) / 2, (dims + thickness) / 2)
    outline = fast_boolean(outline_outer, outline_inner, "not", layer=layer)
    return outline


# def chip(design_func, dims=DEFAULT_DIMS, name=None):
#     cell = Cell('Chip_{}'.format(name) if name else 'Chip')

# TODO: add snake gap using multiple snake calls, put gap between them on chip


@click.group()
def cli():
    pass


if __name__ == "__main__":
    cli()
