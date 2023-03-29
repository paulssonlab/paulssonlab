#!/usr/bin/env python
import numpy as np
import gdstk
from gdstk import Reference, rectangle, boolean, Curve, Polygon
from paulssonlab.microfluidics_design.geometry import (
    Cell,
    ellipse,
    cross,
    qr_target,
    mirror,
    mirror_refs,
    align,
    flatten_or_merge,
)
from paulssonlab.microfluidics_design.text import text as _text
from paulssonlab.microfluidics_design.util import make_odd, memoize, get_uuid, write_gds
import paulssonlab.microfluidics_design.hamming as hamming
import bitarray.util
from functools import partial
from cytoolz import compose
from itertools import product
import numbers

CURVE_TOLERANCE = 0.2
REFERENCE_LAYER = 0
TRENCH_LAYER = 1
FEEDING_CHANNEL_LAYER = 2
DEFAULT_DIMS = np.array([23e3, 13e3])


def text(
    s,
    size,
    position=(0, 0),
    horizontal_alignment=None,
    vertical_alignment=None,
    prerotate_alignment=None,
    angle=0,
    **kwargs,
):
    objs = _text(s, size, position=(0, 0), **kwargs)
    if prerotate_alignment is not None:
        objs = align(objs, position=(0, 0), horizontal=prerotate_alignment)
    if angle == 0:
        pass
    elif angle == -np.pi / 2:
        for ref in objs:
            ref.rotation += np.deg2rad(90)
            ref.origin = np.array(ref.origin)[::-1] * np.array([-1, 1])
    elif angle == np.pi / 2:
        for ref in objs:
            ref.rotation -= np.deg2rad(90)
            ref.origin = np.array(ref.origin)[::-1] * np.array([1, -1])
    else:
        raise NotImplementedError
    objs = align(
        objs,
        position=position,
        horizontal=horizontal_alignment,
        vertical=vertical_alignment,
    )
    return objs


def _compute_lane_split(split, max_lanes, gap_lanes=0):
    """Given a split specification and maximum number of lanes, generates a
    list of lanes per snake.

    Parameters
    ----------
    split : Union[int, Iterable[int]], optional
        If an integer, specifies the number of snakes of approximately-equal size. If a
        tuple, specifies the number of lanes for each snake.
    max_lanes : int
        Maximum number of lanes that can fit on the chip.
    gap_lanes : int, optional
        Number of lanes' worth of blank space to leave between each snake.

    Returns
    -------
    split
        An array specifying how many lanes each snake contains.

    Raises
    ------
    ValueError
        Raised if 'split' is an invalid split specification.
    """
    if split is None:
        split = 1
    if isinstance(split, numbers.Integral):
        max_lanes -= gap_lanes * (split - 1)
        good_split = False
        for lanes_per_snake in (int(np.around(max_lanes / split)), max_lanes // split):
            lanes_per_snake = make_odd(lanes_per_snake)
            remainder_lanes_per_snake = make_odd(
                max_lanes - (split - 1) * lanes_per_snake
            )
            new_split = (lanes_per_snake,) * (split - 1) + (remainder_lanes_per_snake,)
            # TODO: should this be s >= 1?
            if all(s > 1 for s in new_split):
                good_split = True
                split = np.array(new_split)
                break
        if not good_split:
            raise ValueError("bad split: {split}".format(",".join(split)))
    else:
        split = np.asarray(split)
        max_lanes -= gap_lanes * (len(split) - 1)
        if isinstance(split[0], numbers.Integral):
            if np.sum(split) > max_lanes:
                raise ValueError(
                    "total lanes desired {} is greater than maximum number of lanes {}".format(
                        sum(split), max_lanes
                    )
                )
            if np.any(split % 2 == 0):
                raise ValueError("number of lanes per snake must be odd")
        else:
            raise ValueError(
                "snake splitting spec must be a integer or sequence of integers"
            )
    return split


def manifold_snake(
    dims=DEFAULT_DIMS,
    manifold_split=1,
    snake_split=None,
    add_remainder=False,
    lanes_per_snake=1,
    manifold_width=200,
    manifold_input_margin=2e3,
    manifold_bend_margin=0.2e3,
    manifold_margin=200,
    manifold_bend_radius=200,
    manifold_round_radius=True,
    manifold_input_style="bend-out",
    border_margin=600,
    port_margin=600,
    trench_width=1.5,
    trench_length=35,
    trench_fc_overlap=None,
    trench_margin=0.5e3,
    trench_gap=20,
    trench_spacing=2,
    manifold_trench_params=None,
    feeding_channel_width=90,
    port_radius=200,
    port=False,
    port_wayfinder=True,
    port_wayfinder_margin=400,
    port_wayfinder_length=200,
    port_wayfinder_width=100,
    port_wayfinder_orientations=("left", "top", "bottom"),
    registration_marks=False,
    registration_mark_barcodes=False,
    barcode_num_bits=11,
    barcode_rows=4,
    barcode_columns=4,
    mark_size=1,
    mark_spacing=1,
    chip_id=None,
    ticks=False,
    tick_length=5,
    tick_margin=5,
    tick_period=25,
    tick_font_size=None,
    tick_labels=False,
    trenches=True,
    flatten_feeding_channel=False,
    merge_feeding_channel=True,
    feeding_channel_layer=FEEDING_CHANNEL_LAYER,
    trench_layer=TRENCH_LAYER,
    name=None,
    label_font_size=600,
    label_margin=100,
    label=True,
):
    if manifold_input_style not in ("u-turn", "bend-out", "bend-in"):
        raise NotImplementedError
    if name is None:
        name = get_uuid()
    if trench_fc_overlap is None:
        trench_fc_overlap = min(trench_length, feeding_channel_width / 3)
    if tick_font_size is None:
        tick_font_size = tick_length * 2
    if manifold_round_radius is True:
        manifold_round_radius = feeding_channel_width
    if manifold_bend_margin is True:
        manifold_bend_margin = manifold_round_radius
    if manifold_round_radius > manifold_bend_margin:
        # TODO: eliminate this constraint by improving implementation of manifold_round_radius
        raise ValueError(
            "manifold_bend_margin needs to be at least manifold_round_radius"
        )
    # define some dimensions
    effective_trench_length = trench_length + trench_gap / 2
    inner_snake_bend_radius = effective_trench_length
    outer_snake_bend_radius = feeding_channel_width + inner_snake_bend_radius
    if manifold_input_style == "u-turn":
        # ensure port has enough margin on the inner side
        manifold_bend_radius = max(
            manifold_bend_radius, (port_radius + port_margin) / 2
        )
        top_margin = (
            border_margin + manifold_width + manifold_bend_radius + manifold_bend_margin
        )
        horizontal_margin = (
            port_margin
            + port_radius
            + manifold_width / 2
            + 2 * manifold_bend_radius
            + manifold_width
            + manifold_margin
        )
    elif manifold_input_style in ("bend-out", "bend-in"):
        if manifold_input_style == "bend-out":
            horizontal_margin = (
                port_margin
                + 2 * port_radius
                + manifold_input_margin
                + manifold_bend_radius
                + manifold_width
                + manifold_margin
            )
        elif manifold_input_style == "bend-in":
            # ensure port has enough margin on inner side
            manifold_bend_margin = max(
                manifold_bend_margin,
                port_radius + port_margin - manifold_bend_radius - manifold_width / 2,
            )
            horizontal_margin = border_margin + manifold_width + manifold_margin
        top_margin = (
            max(port_margin + port_radius, border_margin + manifold_width / 2)
            + manifold_width / 2
            + manifold_bend_radius
            + manifold_bend_margin
        )
    bottom_margin = top_margin
    if label:
        # only make room for label text if necessary
        total_label_margin = 2 * label_margin + label_font_size
        top_margin = max(top_margin, total_label_margin)
    lane_fc_dims = np.array(
        [
            dims[0] - 2 * horizontal_margin - 2 * outer_snake_bend_radius,
            feeding_channel_width,
        ]
    )
    # split
    lane_height = feeding_channel_width + 2 * effective_trench_length
    max_lanes = int((dims[1] - top_margin - bottom_margin) // lane_height)
    if lanes_per_snake:
        if lanes_per_snake <= 0 or lanes_per_snake % 2 == 0:
            raise ValueError("lanes_per_snake must be positive and odd")
        if snake_split is not None:
            raise ValueError("cannot specify both snake_split and lanes_per_snake")
        snake_split = (lanes_per_snake,) * int(max_lanes // lanes_per_snake)
        if add_remainder:
            remainder = max_lanes % lanes_per_snake
            if remainder % 2 == 0:
                remainder -= 1
            if remainder > 0:
                snake_split += (remainder,)
    snake_split = _compute_lane_split(snake_split, max_lanes)
    num_snakes = len(snake_split)
    snake_split_cum = np.concatenate(((0,), np.cumsum(snake_split)))
    if np.isscalar(manifold_split):
        num_manifolds = manifold_split
        lanes_per_input = int(num_snakes // num_manifolds)
        manifold_split = (lanes_per_input,) * num_manifolds
        remainder = num_snakes % num_manifolds
        manifold_split = np.array(manifold_split)
        if remainder > 0:
            # evenly distribute remainder starting from top input
            manifold_split[:remainder] += 1
    else:
        if np.sum(manifold_split) != num_snakes:
            raise ValueError(
                f"manifold_split must sum to {num_snakes} (total number of snakes)"
            )
    num_manifolds = len(manifold_split)
    manifold_split_cum = np.concatenate(((0,), np.cumsum(manifold_split)))
    num_lanes = np.sum(snake_split)
    lanes_per_input = np.array(
        [
            snake_split[
                manifold_split_cum[:-1][idx] : manifold_split_cum[1:][idx]
            ].sum()
            for idx in range(len(manifold_split_cum) - 1)
        ]
    )
    # allow different sets of trench parameters
    if manifold_trench_params is not None:
        if len(manifold_trench_params) != len(manifold_split):
            raise ValueError(
                f"manifold_trench_params should have the same length as the number of manifolds ({num_manifolds})"
            )
    else:
        manifold_trench_params = [dict() for _ in range(num_manifolds)]
    # root cell
    snake_cell = Cell(f"Snake-{name}")
    # label text
    label_position = (0, dims[1] / 2 - label_margin - label_font_size)
    if label:
        snake_cell.add(
            *text(
                name,
                label_font_size,
                position=label_position,
                horizontal_alignment="center",
                layer=feeding_channel_layer,
            )
        )
    # feeding channel
    snake_fc_cell, lane_ys = _snake_feeding_channel(
        name=name,
        split=snake_split,
        lane_fc_dims=lane_fc_dims,
        effective_trench_length=effective_trench_length,
        port_offset=outer_snake_bend_radius,
        port_radius=port_radius,
        port=False,
        port_wayfinder=False,
        port_wayfinder_margin=None,
        port_wayfinder_length=None,
        port_wayfinder_width=None,
        layer=feeding_channel_layer,
    )
    # manifolds
    for idx in range(len(manifold_split_cum) - 1):
        left_port_lanes = snake_split_cum[:-1][
            manifold_split_cum[idx] : manifold_split_cum[idx + 1]
        ]
        right_port_lanes = (snake_split_cum[1:] - 1)[
            manifold_split_cum[idx] : manifold_split_cum[idx + 1]
        ]
        snake_manifold_cell = _manifold(
            name=name,
            dims=dims,
            lane_fc_dims=lane_fc_dims,
            effective_trench_length=effective_trench_length,
            lane_ys=lane_ys,
            left_port_lanes=left_port_lanes,
            right_port_lanes=right_port_lanes,
            manifold_split_cum=manifold_split_cum,
            feeding_channel_width=feeding_channel_width,
            manifold_width=manifold_width,
            manifold_input_margin=manifold_input_margin,
            manifold_bend_margin=manifold_bend_margin,
            manifold_bend_radius=manifold_bend_radius,
            manifold_round_radius=manifold_round_radius,
            manifold_input_style=manifold_input_style,
            port_margin=port_margin,
            port_radius=port_radius,
            port=port,
            port_wayfinder=port_wayfinder,
            port_wayfinder_margin=port_wayfinder_margin,
            port_wayfinder_length=port_wayfinder_length,
            port_wayfinder_width=port_wayfinder_width,
            port_wayfinder_orientations=port_wayfinder_orientations,
            feeding_channel_layer=feeding_channel_layer,
        )
        snake_fc_cell.add(Reference(snake_manifold_cell, (0, 0)))
    flatten_or_merge(
        snake_fc_cell,
        flatten=flatten_feeding_channel,
        merge=merge_feeding_channel,
        layer=feeding_channel_layer,
    )
    y_offset = (bottom_margin - top_margin) / 2
    snake_cell.add(Reference(snake_fc_cell, (0, y_offset)))
    # add y_offset to lane_ys so they are absolute, all remaining geometry shouldn't
    # need to correct for y_offset
    lane_ys += y_offset
    # trenches
    if trenches:
        trench_active_width = lane_fc_dims[0] - trench_margin - trench_width
        trench_xs = []
        for idx in range(num_manifolds):
            trench_name = f"{name}-T{idx}"
            trench_lane_ys = lane_ys[
                snake_split_cum[manifold_split_cum[idx]] : snake_split_cum[
                    manifold_split_cum[idx + 1]
                ]
            ]
            snake_trenches_cell, manifold_trench_xs = _snake_trenches(
                **{
                    **dict(
                        trench_active_width=trench_active_width,
                        trench_width=trench_width,
                        trench_spacing=trench_spacing,
                        trench_length=trench_length,
                        trench_gap=trench_gap,
                        trench_fc_overlap=trench_fc_overlap,
                        feeding_channel_width=feeding_channel_width,
                        registration_marks=registration_marks,
                        registration_mark_barcodes=registration_mark_barcodes,
                        barcode_num_bits=barcode_num_bits,
                        barcode_rows=barcode_rows,
                        barcode_columns=barcode_columns,
                        mark_size=mark_size,
                        mark_spacing=mark_spacing,
                        chip_id=chip_id,
                        ticks=ticks,
                        tick_labels=tick_labels,
                        tick_margin=tick_margin,
                        tick_length=tick_length,
                        tick_period=tick_period,
                        tick_font_size=tick_font_size,
                        lane_ys=trench_lane_ys,
                        name=trench_name,
                        layer=trench_layer,
                    ),
                    **manifold_trench_params[idx],
                }
            )
            trench_xs.append(manifold_trench_xs)
            snake_cell.add(Reference(snake_trenches_cell, (0, 0)))
    trenches_per_set = np.array([len(xs) for xs in trench_xs])
    trenches_per_input = trenches_per_set * 2 * lanes_per_input
    num_trenches = trenches_per_input.sum()
    lane_length = lane_fc_dims[0]
    metadata = {
        k: v
        for k, v in locals().items()
        if k
        in (
            "num_lanes",
            "trenches_per_set",
            "lanes_per_input",
            "num_trenches",
            "trenches_per_input",
            "feeding_channel_width",
            "trench_gap",
            "trench_length",
            "manifold_width",
            "outer_snake_bend_radius",
            "lane_ys",
            "snake_split",
            "snake_split_cum",
            "manifold_split",
            "manifold_split_cum",
            "trench_xs",
            "lane_length",
        )
    }
    # TODO:
    lane_with_trenches_length = np.array(
        [xs[-1] - xs[0] + trench_width for xs in trench_xs]
    )
    fov_origin_x = np.array([xs[0] - trench_width / 2 for xs in trench_xs])
    selected_trench_region = np.argmax(lane_with_trenches_length)
    metadata["lane_with_trenches_length"] = lane_with_trenches_length[
        selected_trench_region
    ]
    metadata["fov_origin_x"] = fov_origin_x[selected_trench_region]
    metadata["fov_origin_y"] = lane_ys[0] + feeding_channel_width / 2 + trench_length
    return snake_cell, metadata


def snake(
    dims=DEFAULT_DIMS,
    split=1,
    border_margin=600,
    port_margin=600,
    port_input_margin=600,
    trench_width=1.5,
    trench_length=35,
    trench_fc_overlap=None,
    trench_margin=0.5e3,
    trench_gap=20,
    gap_lanes=0,
    trench_spacing=2,
    feeding_channel_width=90,
    port_radius=200,
    port=False,
    port_wayfinder=True,
    port_wayfinder_margin=400,
    port_wayfinder_length=200,
    port_wayfinder_width=100,
    registration_marks=False,
    registration_mark_barcodes=False,
    barcode_num_bits=11,
    barcode_rows=4,
    barcode_columns=4,
    mark_size=1,
    mark_spacing=1,
    chip_id=None,
    ticks=False,
    tick_length=5,
    tick_margin=5,
    tick_period=25,
    tick_font_size=None,
    tick_labels=False,
    trenches=True,
    flatten_feeding_channel=False,
    merge_feeding_channel=True,
    feeding_channel_layer=FEEDING_CHANNEL_LAYER,
    trench_layer=TRENCH_LAYER,
    label_font_size=600,
    label_margin=100,
    label=True,
    name=None,
):
    """Summary.

    Parameters
    ----------
    dims : Tuple[float, float], optional
        A tuple of `(height, width)` specifying the dimensions of the chip footprint (in
        microns).
    split : Union[int, Iterable[int]], optional
        If an integer, specifies the number of snakes of approximately-equal size. If a
        tuple, specifies the number of lanes for each snake.
    top_margin : float, optional
        Empty space (in microns) between the top dead end of the topmost set of trenches
        and the bottom edge of the chip.
    bottom_margin : float, optional
        Empty space (in microns) between the bottom dead end of the bottommost set of
        trenches and the bottom edge of the chip.
    port_margin : float, optional
        Empty space (in microns) between the outside of each port and the left/right
        edges of the chip.
    trench_width : float, optional
        Width (in microns) of the trenches.
    trench_length : float, optional
        Length (in microns) of the trenches.
    trench_fc_overlap : float, optional
        Length (in microns) that trenches should extend into the feeding channel. This
        allows the design to tolerate misalignment between the feeding channel and
        trench layers.
    trench_margin : float, optional
        Trenches will be excluded from a region this length (in microns) on each end of
        each straight section of feeding channel.
    trench_gap : float, optional
        Distance (in microns) between the dead end of one set of trenches and the dead
        ends of the next set of trenches.
    gap_lanes : int, optional
        If nonzero, the number of lanes that are left without a feeding channel. Setting
        this to 1 could decrease the probability of flow between adjacent snakes due to
        defects bridging two snakes.
    trench_spacing : float, optional
        Empty space (in microns) between the edge of one trench and the edge of the next
        trench.
    feeding_channel_width : float, optional
        Width (in microns) of feeding channel.
    port_radius : float, optional
        Radius (in microns) of the ports.
    ticks : bool, optional
        Whether to include small rectangular tick marks above the dead-ends of the
        trenches every `tick_period` trenches.
    tick_length : float, optional
        Length (in microns) of the tick marks.
    tick_margin : float, optional
        Distance (in microns) between the dead end of a trench and the tick marks.
    tick_period : int, optional
        Number of trenches between tick marks.
    tick_font_size : None, optional
        Font size (in microns) of the trench tick numbering.
    tick_labels : bool, optional
        Whether to label every tick mark with the trench number.
    trenches : bool, optional
        Whether to output trenches or not.
    flatten_feeding_channel : bool, optional
        If False, each bend or straight section of a feeding channel will be specified
        as a cell reference. If True, each polygon of the feeding channel will be output
        individually.
    merge_feeding_channel : bool, optional
        If False, the feeding channel bends and the straight rectangular sections will
        be output as separate polygons. If True, each snake's feeding channel will be
        output as a single polygon. The MLA150 software seems to prefer merged feeding
        channels.
    feeding_channel_layer : int, optional
        Feeding channel layer number.
    trench_layer : int, optional
        Trench layer number.
    label_font_size : float, optional
        Font size (in microns) of the label text.
    label_margin : float, optional
        Distance between the top edge of chip and the top of the label text.
    label : bool, optional
        If True, a label is written along the top edge of the chip.
    name : str, optional
        Unique identifier used for cell names. If None is given, an alphanumeric UUID
        will be generated.

    Returned
    ------------------
    snake_cell
        GDS cell containing the chip design.
    metadata : dict
        Dict containing the design parameters.
    """
    if name is None:
        name = get_uuid()
    if trench_fc_overlap is None:
        trench_fc_overlap = min(trench_length, feeding_channel_width / 3)
    if tick_font_size is None:
        tick_font_size = tick_length * 2
    effective_trench_length = trench_length + trench_gap / 2
    inner_snake_bend_radius = effective_trench_length
    outer_snake_bend_radius = feeding_channel_width + inner_snake_bend_radius
    top_margin = max(
        border_margin, port_margin + port_radius - feeding_channel_width / 2
    )
    bottom_margin = top_margin
    if label:
        # only make room for label text if necessary
        total_label_margin = 2 * label_margin + label_font_size
        top_margin = max(top_margin, total_label_margin)
    horizontal_margin = port_margin + 2 * port_radius + port_input_margin
    lane_fc_dims = np.array(
        [
            dims[0] - 2 * horizontal_margin - outer_snake_bend_radius,
            feeding_channel_width,
        ]
    )
    lane_height = feeding_channel_width + 2 * effective_trench_length
    max_lanes = int((dims[1] - top_margin - bottom_margin) // lane_height)
    split = _compute_lane_split(split, max_lanes, gap_lanes=gap_lanes)
    num_lanes = np.sum(split)
    lanes_per_input = np.array([num_lanes])
    trenches_per_input = np.array([num_trenches])
    snake_cell = Cell(f"Snake-{name}")
    # label text
    label_position = (0, dims[1] / 2 - label_margin - label_font_size)
    if label:
        snake_cell.add(
            *text(
                name,
                label_font_size,
                position=label_position,
                horizontal_alignment="center",
                layer=feeding_channel_layer,
            )
        )
    port_offset = dims[0] / 2 - lane_fc_dims[0] / 2 - port_radius - port_margin
    snake_fc_cell, lane_ys = _snake_feeding_channel(
        name=name,
        split=split,
        lane_fc_dims=lane_fc_dims,
        effective_trench_length=effective_trench_length,
        port_offset=port_offset,
        port_radius=port_radius,
        port=port,
        port_wayfinder=port_wayfinder,
        port_wayfinder_margin=port_wayfinder_margin,
        port_wayfinder_length=port_wayfinder_length,
        port_wayfinder_width=port_wayfinder_width,
        gap_lanes=gap_lanes,
        layer=feeding_channel_layer,
        flatten_feeding_channel=flatten_feeding_channel,
        merge_feeding_channel=merge_feeding_channel,
    )
    flatten_or_merge(
        snake_fc_cell,
        flatten=flatten_feeding_channel,
        merge=merge_feeding_channel,
        layer=feeding_channel_layer,
    )
    y_offset = (bottom_margin - top_margin) / 2
    snake_cell.add(Reference(snake_fc_cell, (0, y_offset)))
    if trenches:
        trench_active_width = lane_fc_dims[0] - trench_margin - trench_width
        snake_trenches_cell, trench_xs = _snake_trenches(
            trench_active_width=trench_active_width,
            trench_width=trench_width,
            trench_spacing=trench_spacing,
            trench_length=trench_length,
            trench_gap=trench_gap,
            trench_fc_overlap=trench_fc_overlap,
            feeding_channel_width=feeding_channel_width,
            registration_marks=registration_marks,
            registration_mark_barcodes=registration_mark_barcodes,
            barcode_num_bits=barcode_num_bits,
            barcode_rows=barcode_rows,
            barcode_columns=barcode_columns,
            mark_size=mark_size,
            mark_spacing=mark_spacing,
            chip_id=chip_id,
            ticks=ticks,
            tick_labels=tick_labels,
            tick_margin=tick_margin,
            tick_length=tick_length,
            tick_period=tick_period,
            tick_font_size=tick_font_size,
            trench_xs=trench_xs,
            lane_ys=lane_ys,
            name=name,
            layer=trench_layer,
        )
        trenches_per_set = len(trench_xs)
        num_trenches = trenches_per_set * 2 * num_lanes
        snake_cell.add(Reference(snake_trenches_cell, (0, y_offset)))
    lane_length = lane_fc_dims[0]
    metadata = {
        k: v
        for k, v in locals().items()
        if k
        in (
            "num_lanes",
            "trenches_per_set",
            "lanes_per_input",
            "num_trenches",
            "trenches_per_input",
            "split",
            "feeding_channel_width",
            "trench_gap",
            "trench_length",
            "trench_xs",
            "lane_length",
        )
    }
    metadata["lane_with_trenches_length"] = trench_xs[-1] - trench_xs[0] + trench_width
    metadata["fov_origin_x"] = trench_xs[0] - trench_width / 2
    metadata["fov_origin_y"] = lane_ys[0] + feeding_channel_width / 2 + trench_length
    snake_length = split * lane_length
    metadata["snake_length"] = snake_length
    return snake_cell, metadata


def _snake_feeding_channel(
    name,
    split,
    lane_fc_dims,
    effective_trench_length,
    port_offset,
    port_radius,
    port,
    port_wayfinder,
    port_wayfinder_margin,
    port_wayfinder_length,
    port_wayfinder_width,
    gap_lanes=0,
    layer=FEEDING_CHANNEL_LAYER,
    flatten_feeding_channel=False,
    merge_feeding_channel=True,
):
    feeding_channel_width = lane_fc_dims[1]
    lane_height = feeding_channel_width + 2 * effective_trench_length
    inner_snake_bend_radius = effective_trench_length
    outer_snake_bend_radius = feeding_channel_width + inner_snake_bend_radius
    lane_cell = Cell(f"Lane-{name}")
    lane_fc = rectangle(-lane_fc_dims / 2, lane_fc_dims / 2, layer=layer)
    lane_cell.add(lane_fc)
    bend_cell = Cell(f"Feeding Channel Bend-{name}")
    bend = ellipse(
        (0, 0),
        outer_snake_bend_radius,
        inner_radius=inner_snake_bend_radius,
        initial_angle=3 / 2 * np.pi,
        final_angle=5 / 2 * np.pi,
        layer=layer,
    )
    bend_cell.add(bend)
    port_cell = Cell(f"Feeding Channel Port-{name}")
    if port and port_radius:
        port = ellipse((port_offset, 0), port_radius, layer=layer)
        port_cell.add(port)
    port_fc = rectangle(
        (0, -feeding_channel_width / 2),
        (port_offset, feeding_channel_width / 2),
        layer=layer,
    )
    port_cell.add(port_fc)
    if port_wayfinder:
        wf = wayfinder(
            radius=port_radius + port_wayfinder_margin,
            length=port_wayfinder_length,
            width=port_wayfinder_width,
            orientations=("left", "top", "bottom"),
        )
        port_cell.add(Reference(wf, (port_offset, 0)))
    max_lanes = sum(split) + gap_lanes * (len(split) - 1)
    last_lane_y = ((max_lanes - 1) * lane_height) / 2
    lane_ys = np.linspace(last_lane_y, -last_lane_y, max_lanes)
    lane_mask = np.full(len(lane_ys), True)
    skipped_lanes = np.cumsum(split + gap_lanes) - gap_lanes
    for offset in range(gap_lanes):
        lane_mask[skipped_lanes[:-1] + offset] = False
    lane_mask[skipped_lanes[-1] :] = False
    lane_ys = lane_ys[lane_mask]
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
    snake_fc_cell = Cell(f"Snake Feeding Channel-{name}")
    for y in lane_ys:
        snake_fc_cell.add(Reference(lane_cell, (0, y)))
    for lane in right_bend_lanes:
        snake_fc_cell.add(
            Reference(bend_cell, (lane_fc_dims[0] / 2, lane_ys[lane] - lane_height / 2))
        )
    for lane in left_bend_lanes:
        snake_fc_cell.add(
            Reference(
                bend_cell,
                (-lane_fc_dims[0] / 2, lane_ys[lane] - lane_height / 2),
                rotation=np.deg2rad(180),
            )
        )
    for lane in right_port_lanes:
        snake_fc_cell.add(Reference(port_cell, (lane_fc_dims[0] / 2, lane_ys[lane])))
    for lane in left_port_lanes:
        snake_fc_cell.add(
            Reference(
                port_cell,
                (-lane_fc_dims[0] / 2, lane_ys[lane]),
                rotation=np.deg2rad(180),
            )
        )
    return snake_fc_cell, lane_ys


def _manifold(
    name,
    dims,
    lane_fc_dims,
    effective_trench_length,
    lane_ys,
    left_port_lanes,
    right_port_lanes,
    manifold_split_cum,
    feeding_channel_width,
    manifold_width,
    manifold_input_margin,
    manifold_bend_margin,
    manifold_bend_radius,
    manifold_round_radius,
    manifold_input_style,
    port_margin,
    port_radius,
    port,
    port_wayfinder,
    port_wayfinder_margin,
    port_wayfinder_length,
    port_wayfinder_width,
    port_wayfinder_orientations,
    feeding_channel_layer,
):
    inner_snake_bend_radius = effective_trench_length
    outer_snake_bend_radius = feeding_channel_width + inner_snake_bend_radius
    if manifold_round_radius:
        rounded_corner = Cell(f"Snake-round-{name}")
        rounded_curve = Curve((0, 0), tolerance=CURVE_TOLERANCE)
        rounded_curve.vertical(manifold_round_radius)
        rounded_curve.arc(manifold_round_radius, 0, -1 / 2 * np.pi, 0)
        rounded_corner.add(Polygon(rounded_curve.points(), layer=feeding_channel_layer))
    manifold_cell = Cell(f"Snake-Manifold-{name}")
    for flip in (1, -1):
        port_lanes = left_port_lanes if flip == -1 else right_port_lanes
        manifold_lane_ys = lane_ys[port_lanes][::-flip]
        manifold_input_bend_y = manifold_lane_ys[0] - flip * (
            feeding_channel_width / 2 + manifold_bend_margin
        )
        if manifold_input_style == "u-turn":
            wf_base_rotation = np.pi
            port_x = -dims[0] / 2 + port_margin + port_radius
            manifold_input_bend_x = port_x + manifold_width / 2 + manifold_bend_radius
            manifold_left_x = manifold_input_bend_x + manifold_bend_radius
            port_y = manifold_input_bend_y + flip * (
                manifold_input_margin + port_radius
            )
            manifold_bend_angles = (0, -flip * np.pi)
        elif manifold_input_style == "bend-out":
            wf_base_rotation = np.pi / 2
            port_x = -dims[0] / 2 + port_margin + port_radius
            manifold_input_bend_x = port_x + port_radius + manifold_input_margin
            manifold_left_x = manifold_input_bend_x + manifold_bend_radius
            port_y = manifold_lane_ys[0] - flip * (
                feeding_channel_width / 2
                + manifold_bend_margin
                + manifold_bend_radius
                + manifold_width / 2
            )
            manifold_bend_angles = (flip + np.array([1, 2])) / 2 * np.pi
        elif manifold_input_style == "bend-in":
            wf_base_rotation = np.pi / 2
            manifold_left_x = -dims[0] / 2 + border_margin
            manifold_input_bend_x = (
                manifold_left_x + manifold_width + manifold_bend_radius
            )
            port_x = manifold_input_bend_x + manifold_input_margin + port_radius
            # TODO: + top_margin?
            port_y = manifold_lane_ys[0] - flip * (
                feeding_channel_width / 2
                + manifold_bend_margin
                + manifold_bend_radius
                + manifold_width / 2
            )
            manifold_bend_angles = (flip + np.array([3, 2])) / 2 * np.pi
        manifold_top_y = manifold_lane_ys[0] - flip * feeding_channel_width / 2
        manifold_taper_y = manifold_lane_ys[-1] + flip * (
            feeding_channel_width / 2 - manifold_width
        )
        if port:
            snake_manifold_cell.add(
                ellipse(
                    (-flip * port_x, port_y),
                    port_radius,
                    layer=feeding_channel_layer,
                )
            )
        if port_wayfinder:
            wf = wayfinder(
                radius=port_radius + port_wayfinder_margin,
                length=port_wayfinder_length,
                width=port_wayfinder_width,
                orientations=port_wayfinder_orientations,
            )
            manifold_cell.add(
                Reference(
                    wf,
                    (-flip * port_x, port_y),
                    rotation=wf_base_rotation + np.deg2rad(90 * flip),
                )
            )
        manifold_cell.add(
            ellipse(
                (-flip * manifold_input_bend_x, manifold_input_bend_y),
                manifold_bend_radius + manifold_width,
                inner_radius=manifold_bend_radius,
                initial_angle=manifold_bend_angles[0],
                final_angle=manifold_bend_angles[1],
                layer=feeding_channel_layer,
            )
        )
        if manifold_input_style in ("bend-out", "bend-in"):
            manifold_cell.add(
                rectangle(
                    (-flip * port_x, port_y + manifold_width / 2),
                    (-flip * manifold_input_bend_x, port_y - manifold_width / 2),
                    layer=feeding_channel_layer,
                )
            )
        elif manifold_input_style == "u-turn":
            manifold_cell.add(
                rectangle(
                    (-flip * (port_x - manifold_width / 2), manifold_input_bend_y),
                    (-flip * (port_x + manifold_width / 2), port_y),
                    layer=feeding_channel_layer,
                )
            )
        # manifold bend margin
        manifold_cell.add(
            rectangle(
                (-flip * manifold_left_x, manifold_input_bend_y),
                (-flip * (manifold_left_x + manifold_width), manifold_top_y),
                layer=feeding_channel_layer,
            )
        )
        # manifold
        manifold_cell.add(
            rectangle(
                (-flip * manifold_left_x, manifold_top_y),
                (-flip * (manifold_left_x + manifold_width), manifold_taper_y),
                layer=feeding_channel_layer,
            )
        )
        if manifold_round_radius:
            for y in manifold_lane_ys[:-1]:
                manifold_cell.add(
                    Reference(
                        rounded_corner,
                        (
                            -flip * (manifold_left_x + manifold_width),
                            y + flip * feeding_channel_width / 2,
                        ),
                        rotation=np.deg2rad(90 * (3 + flip)),
                    )
                )
            for y in manifold_lane_ys:
                manifold_cell.add(
                    Reference(
                        rounded_corner,
                        (
                            -flip * (manifold_left_x + manifold_width),
                            y - flip * feeding_channel_width / 2,
                        ),
                        rotation=np.deg2rad(90 * (4 + flip)),
                    )
                )
        for y in manifold_lane_ys:
            manifold_cell.add(
                rectangle(
                    (
                        -flip * (manifold_left_x + manifold_width),
                        y + feeding_channel_width / 2,
                    ),
                    (
                        -flip * (-lane_fc_dims[0] / 2 - outer_snake_bend_radius),
                        y - feeding_channel_width / 2,
                    ),
                    layer=feeding_channel_layer,
                )
            )
        manifold_taper_angle = np.pi * (1 / 2 - flip)
        manifold_cell.add(
            ellipse(
                (-flip * (manifold_left_x + manifold_width), manifold_taper_y),
                manifold_width,
                initial_angle=manifold_taper_angle,
                final_angle=manifold_taper_angle + flip * np.pi,
                layer=feeding_channel_layer,
            )
        )
    return manifold_cell


@memoize
def _barcode(
    ary,
    mark_size,
    mark_spacing,
    zero_symbol=None,
    one_symbol=None,
    columns=None,
    rows=None,
    layer=None,
):
    if zero_symbol is None:
        zero_symbol = Polygon([(1 / 2, -1 / 2), (1 / 2, 1 / 2), (-1 / 2, -1 / 2)])
    elif zero_symbol is not False:
        zero_symbol = gdstk.copy(zero_symbol)
    if one_symbol is None:
        one_symbol = rectangle((-1 / 2, -1 / 2), (1 / 2, 1 / 2))
    elif one_symbol is not False:
        one_symbol = gdstk.copy(one_symbol)
    for symbol in (zero_symbol, one_symbol):
        if symbol is False:
            continue
        symbol.scale(mark_size)
        symbol.layers = [layer for l in symbol.layers]
    if columns is None and rows is None:
        raise ValueError("either columns or rows must be specified")
    N = len(ary)
    if columns is None:
        rows = int(np.ceil(N / columns))
    elif rows is None:
        columns = int(np.ceil(N / rows))
    padding = rows * columns - N
    ary2 = bitarray.util.zeros(padding)
    ary2.extend(ary)
    number = bitarray.util.ba2int(ary)  # TODO: does not handle different endiannesses
    cell = Cell(f"Barcode-m{mark_size}-s{mark_spacing}-{padding}-{number}")
    mark_pitch = mark_size + mark_spacing
    y_offset = -((rows - 1) * mark_pitch) / 2
    for column in range(columns):
        for row in range(rows):
            x = mark_pitch * column + mark_size / 2
            y = y_offset + mark_pitch * row
            idx = row + rows * column
            if ary2[idx]:
                if one_symbol is not False:
                    cell.add(gdstk.copy(one_symbol, x, y))
            else:
                if zero_symbol is not False:
                    cell.add(gdstk.copy(zero_symbol, x, y))
    return cell


def _snake_trenches(
    trench_active_width,
    trench_width,
    trench_spacing,
    trench_length,
    trench_gap,
    trench_fc_overlap,
    feeding_channel_width,
    registration_marks,
    registration_mark_barcodes,
    mark_size,
    mark_spacing,
    barcode_num_bits,
    barcode_rows,
    barcode_columns,
    chip_id,
    ticks,
    tick_labels,
    tick_margin,
    tick_length,
    tick_period,
    tick_font_size,
    lane_ys,
    name,
    layer=TRENCH_LAYER,
):
    if ticks and (registration_marks or registration_mark_barcodes):
        raise ValueError("cannot draw both ticks and registration marks")
    trench_xs = np.arange(
        -trench_active_width / 2,
        trench_active_width / 2,
        trench_width + trench_spacing,
    )
    trenches_per_set = len(trench_xs)
    lane_gap_offset_y = feeding_channel_width / 2 + trench_length + trench_gap / 2
    mark_pitch = mark_size + mark_spacing
    column_barcode_margin = (
        4 * mark_size + 3 * mark_spacing
    ) / 2 + mark_spacing  # TODO: couple to qr_target arguments, below
    row_barcode_margin = column_barcode_margin + barcode_columns * mark_pitch
    chip_barcode_margin = row_barcode_margin + barcode_columns * mark_pitch
    trenches_per_set = len(trench_xs)
    snake_trenches_cell = Cell(f"Snake Trenches-{name}")
    trench_cell = Cell(f"Trench-{name}")
    trench_cell.add(
        rectangle(
            (-trench_width / 2, -trench_fc_overlap),
            (trench_width / 2, trench_length),
            layer=layer,
        )
    )
    tick_cell = Cell(f"Tick-{name}")
    if ticks:
        tick_cell.add(
            rectangle(
                (-trench_width / 2, tick_margin - trench_gap / 2),
                (trench_width / 2, tick_margin - trench_gap / 2 + tick_length),
                layer=layer,
            )
        )
    elif registration_marks:
        mark_halfwidth = (2 * mark_size + mark_spacing) / 2
        if registration_marks == "qr":
            tick_cell.add(
                *qr_target(
                    mark_size, mark_spacing, 2 * mark_size + mark_spacing, layer=layer
                )
            )
        elif registration == "box" or registration_marks is True:
            tick_cell.add(
                *rectangle(
                    (-mark_halfwidth, -mark_halfwidth),
                    (mark_halfwidth, mark_halfwidth),
                    layer=layer,
                )
            )
        else:
            raise ValueError(f"unknown registration_marks style {registration_marks}")
    tick_xs = trench_xs[::tick_period]
    num_ticks = len(tick_xs)
    if registration_mark_barcodes:
        lane_ys_diff = np.diff(lane_ys)
        uniform_lane_ys = np.all(lane_ys_diff == lane_ys_diff[0])
        if chip_id is not None:
            bits = hamming.encode(
                bitarray.util.int2ba(chip_id, length=barcode_num_bits)
            )
            chip_barcode = _barcode(
                bits,
                mark_size,
                mark_spacing,
                rows=barcode_rows,
                columns=barcode_columns,
                layer=layer,
            )
            if uniform_lane_ys:
                snake_trenches_cell.add(
                    Reference(
                        chip_barcode,
                        (
                            trench_xs[0] + chip_barcode_margin,
                            lane_ys[0] + lane_gap_offset_y,
                        ),
                        column=num_ticks,
                        rows=len(lane_ys),
                        spacing=(
                            tick_period * (trench_width + trench_spacing),
                            lane_ys_diff[0],
                        ),
                    )
                )
            else:
                for y in lane_ys:
                    snake_trenches_cell.add(
                        Reference(
                            chip_barcode,
                            (trench_xs[0] + row_barcode_margin, y + lane_gap_offset_y),
                            columns=num_ticks,
                            rows=1,
                            spacing=(tick_period * (trench_width + trench_spacing), 0),
                        )
                    )
        for tick_idx, x in enumerate(tick_xs):
            bits = hamming.encode(
                bitarray.util.int2ba(tick_idx, length=barcode_num_bits)
            )
            column_barcode = _barcode(
                bits,
                mark_size,
                mark_spacing,
                rows=barcode_rows,
                columns=barcode_columns,
                layer=layer,
            )
            if uniform_lane_ys:
                snake_trenches_cell.add(
                    Reference(
                        column_barcode,
                        (x + column_barcode_margin, lane_ys[0] + lane_gap_offset_y),
                        columns=1,
                        rows=len(lane_ys),
                        spacing=(0, lane_ys_diff[0]),
                    )
                )
            else:
                for y in lane_ys:
                    snake_trenches_cell.add(
                        Reference(
                            column_barcode,
                            (x + column_barcode_margin, y + lane_gap_offset_y),
                        )
                    )
    for lane_idx, y in enumerate(lane_ys):
        snake_trenches_cell.add(
            Reference(
                trench_cell,
                (trench_xs[0], y + feeding_channel_width / 2),
                columns=trenches_per_set,
                rows=1,
                spacing=(trench_width + trench_spacing, 0),
            )
        )
        snake_trenches_cell.add(
            Reference(
                trench_cell,
                (trench_xs[0], y - feeding_channel_width / 2),
                columns=trenches_per_set,
                rows=1,
                spacing=(trench_width + trench_spacing, 0),
                x_reflection=True,
            )
        )
        if ticks or registration_marks or registration_mark_barcodes:
            snake_trenches_cell.add(
                Reference(
                    tick_cell,
                    (trench_xs[0], y + lane_gap_offset_y),
                    columns=num_ticks,
                    rows=1,
                    spacing=(tick_period * (trench_width + trench_spacing), 0),
                )
            )
            if registration_mark_barcodes:
                bits = hamming.encode(
                    bitarray.util.int2ba(lane_idx, length=barcode_num_bits)
                )
                row_barcode = _barcode(
                    bits,
                    mark_size,
                    mark_spacing,
                    rows=barcode_rows,
                    columns=barcode_columns,
                    layer=layer,
                )
                snake_trenches_cell.add(
                    Reference(
                        row_barcode,
                        (trench_xs[0] + row_barcode_margin, y + lane_gap_offset_y),
                        columns=num_ticks,
                        rows=1,
                        spacing=(tick_period * (trench_width + trench_spacing), 0),
                    )
                )
        if tick_labels:
            for tick_idx, x in enumerate(tick_xs):
                tick_idx = (
                    lane_idx * 2 * trenches_per_set + 2 * tick_idx * tick_period + 1
                )
                snake_trenches_cell.add(
                    *text(
                        str(tick_idx),
                        tick_font_size,
                        (
                            x + 2 * trench_width,
                            y + feeding_channel_width / 2 + trench_length + tick_margin,
                        ),
                        layer=layer,
                    )
                )
    return snake_trenches_cell, trench_xs


@memoize
def wayfinder(
    radius=400,
    port=False,
    length=150,
    width=100,
    orientations=None,
    layer=FEEDING_CHANNEL_LAYER,
):
    horizontal_flips = []
    vertical_flips = []
    if orientations is None:
        orientations = ("left", "right", "top", "bottom")
    if "left" in orientations:
        horizontal_flips.append(-1)
    if "right" in orientations:
        horizontal_flips.append(1)
    if "bottom" in orientations:
        vertical_flips.append(-1)
    if "top" in orientations:
        vertical_flips.append(1)
    wayfinder_cell = Cell("Wayfinder")
    for vertical_flip in vertical_flips:
        wayfinder_cell.add(
            rectangle(
                (-width / 2, vertical_flip * radius),
                (width / 2, vertical_flip * (radius + length)),
                layer=layer,
            )
        )
    for horizontal_flip in horizontal_flips:
        wayfinder_cell.add(
            rectangle(
                (horizontal_flip * radius, -width / 2),
                (horizontal_flip * (radius + length), width / 2),
                layer=layer,
            )
        )
    return wayfinder_cell


@memoize
def profilometry_marks(
    dims=np.array([150, 75]),
    columns=3,
    rows=3,
    layers=(FEEDING_CHANNEL_LAYER, TRENCH_LAYER),
    label=True,
):
    """Generates a grid of rectangular profilometry targets.

    Parameters
    ----------
    dims : Tuple[float, float], optional
        A tuple of `(height, width)` specifying the dimensions of the chip footprint (in
        microns).
    columns : int, optional
        Number of columns.
    rows : int, optional
        Number of rows.
    layers : int, optional
        Layer number.
    label : bool, optional
        If True, the layer number is written next to the corresponding targets.

    Returns
    -------
    profilometry_cell
        GDS cell containing the profilometry targets.
    """
    dims = np.array(dims)
    grid_size = np.array([rows, columns])
    grid_dims = dims * grid_size * 2
    profilometry_cell = Cell("Profilometry Marks")
    mark_cells = {layer: Cell(f"Profilometry Mark Layer {layer}") for layer in layers}
    for i, (layer, cell) in enumerate(mark_cells.items()):
        cell.add(rectangle(dims, -dims / 2, layer=layer))
        origin = -grid_dims / 2 + (i - (len(mark_cells) - 1) / 2) * np.array(
            [grid_dims[0], 0]
        )
        profilometry_cell.add(
            Reference(cell, origin, columns=columns, rows=rows, spacing=dims * 2)
        )
        if label:
            profilometry_cell.add(
                *text(
                    str(layer),
                    dims[1],
                    origin - np.array([0, 2 * dims[1]]),
                    layer=layer,
                )
            )
    return profilometry_cell


@memoize
def alignment_cross(length=1e3, thickness=6, layer=TRENCH_LAYER):
    """Makes an alignment cross compatible with the MLA150.

    Parameters
    ----------
    length : float, optional
        Length (in microns) of each appendage of the cross.
    width : float, optional
        Thickness (in microns) of the cross.
    layer : int, optional
        Layer number.

    Returns
    -------
    alignment_cell
        GDS cell containing the alignment cross.
    """
    alignment_cell = Cell("Alignment Cross")
    alignment_cell.add(cross(length, thickness, layer=layer))
    alignment_cell.add(
        rectangle((-3 * length, -length), (-2 * length, length), layer=layer)
    )
    alignment_cell.add(
        rectangle((2 * length, -length), (3 * length, length), layer=layer)
    )
    return alignment_cell


@memoize
def mask_alignment_cross(
    length=300,
    thickness=100,
    cross_spacing=30,
    num_crosses=2,
    bottom_layer=TRENCH_LAYER,
    top_layer=FEEDING_CHANNEL_LAYER,
):
    """Makes an alignment cross compatible with a mask aligner.

    Parameters
    ----------
    length : float, optional
        Length of each individual cross appendage (in microns).
    width : float, optional
        Thickness of each cross (in microns).
    cross_spacing : float, optional
        Spacing (in microns) between crosses.
    num_crosses : int, optional
        Number of crosses.
    bottom_layer : int, optional
        Bottom layer number.
    top_layer : int, optional
        Top layer number.

    Returns
    -------
    alignment_cell
        GDS cell containing the alignment crosses.
    """
    alignment_cell = Cell("Mask Alignment Cross")
    offset_unit = 2 * length + cross_spacing
    box_length = offset_unit * (num_crosses + 1)
    box_corner = np.array([box_length, box_length])
    outer_box = rectangle(-(box_corner + thickness), box_corner + thickness)
    inner_box = rectangle(-box_corner, box_corner)
    box = boolean(outer_box, inner_box, "not", layer=top_layer)
    cross_box_corner = np.array([length, length])
    base_cross_box = rectangle(-cross_box_corner, cross_box_corner, layer=bottom_layer)
    base_cross = cross(length, thickness, layer=top_layer)
    base_cross_box = boolean(base_cross_box, base_cross, "not", layer=bottom_layer)
    alignment_cell.add(base_cross_box)
    alignment_cell.add(base_cross)
    for offset_unit_vector in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        offset_vector = np.array(offset_unit_vector) * offset_unit
        for num_offsets in range(1, num_crosses + 1):
            offset = offset_vector * num_offsets
            alignment_cell.add(gdstk.copy(base_cross_box, *offset))
            alignment_cell.add(gdstk.copy(base_cross, *offset))
    alignment_cell.add(box)
    return alignment_cell


def wafer(
    chips,
    label_right=None,
    label_left=None,
    diameter=76.2e3,
    chip_dims=None,
    chip_margin=1.2e3,
    alignment_mark_position=None,
    alignment_font_size=1000,
    right_font_size=2000,
    left_font_size=1300,
    label=True,
    mask=False,
    feeding_channel_layer=FEEDING_CHANNEL_LAYER,
    trench_layer=TRENCH_LAYER,
    reference_layer=REFERENCE_LAYER,
):
    """Generates a wafer containing the given chips.

    Parameters
    ----------
    chips : List
        List of GDS cells containing chips.
    label_right : str
        Text to write on right side of wafer.
    label_left : str
        Text to write on left side of wafer.
    diameter : float, optional
        Diameter of wafer (in microns).
    chip_dims : Tuple[float, float], optional
        Description
    chip_margin : float, optional
        Description
    alignment_mark_position : float, optional
        Distance (in microns) between both alignment marks.
    alignment_font_size : float, optional
        Size (in microns) of label text for alignment marks.
    right_font_size : float, optional
        Size (in microns) of text.
    left_font_size : float, optional
        Size (in microns) of text.
    label : bool, optional
        If true, label the wafer with `label_right` and `label_left`. If mask is true, additionally label the
        masks with `name` and the layer number outside of the wafer area.
    mask : bool, optional
        If true, mask aligner-compatible alignment crosses are generated. If false,
        MLA150-compatible alignment marks are generated.
    feeding_channel_layer : int, optional
        Feeding channel layer number.
    trench_layer : int, optional
        Trench layer number.
    reference_layer : int, optional
        Layer number for wafer and chip area outline.
    """
    if len(chips) > 6:
        raise Exception("cannot lay out more than six chips on a wafer")
    if chip_dims is None:
        raise Exception("must provide maximum chip dimension")
    main_cell = Cell("main")
    square_corner = np.array([diameter, diameter]) * 1.2 / 2
    square = rectangle(-square_corner, square_corner, layer=reference_layer)
    circle = ellipse((0, 0), diameter / 2)
    wafer_outline = boolean(square, circle, "not", layer=reference_layer)
    main_cell.add(*wafer_outline)
    horizontal_chip_spacing = chip_dims[0] + chip_margin
    vertical_chip_spacing = chip_dims[1] + chip_margin
    vertical_spacings = (-1, 0, 1)
    chips = [*reversed(chips[:3]), *chips[3:]]
    if len(chips) <= 3:
        horizontal_spacings = (0,)
    else:
        horizontal_spacings = (-1, 1)
    for (horizontal, vertical), chip in zip(
        product(horizontal_spacings, vertical_spacings), chips
    ):
        x = horizontal_chip_spacing * horizontal / 2
        y = vertical_chip_spacing * vertical  # * 2 / 3
        main_cell.add(Reference(chip, (x, y)))
    chip_area = np.array(
        (
            len(horizontal_spacings) * chip_dims[0]
            + (len(horizontal_spacings) - 1) * chip_margin,
            len(vertical_spacings) * chip_dims[1]
            + (len(vertical_spacings) - 1) * chip_margin,
        )
    )
    chip_area_corner = chip_area / 2
    profilometry_cell = profilometry_marks(
        layers=(feeding_channel_layer, trench_layer), label=label
    )
    profilometry_bbox = np.array(profilometry_cell.bounding_box())
    profilometry_spacing = np.array(
        [chip_area[0] + 2 * np.abs(profilometry_bbox[:, 0]).max(), 0]
    )
    main_cell.add(
        Reference(
            profilometry_cell,
            np.array(
                (-chip_area[0] / 2 - profilometry_bbox[:, 0].max() - chip_margin, 0)
            ),
        )
    )
    main_cell.add(
        Reference(
            profilometry_cell,
            np.array(
                (chip_area[0] / 2 - profilometry_bbox[:, 0].min() + chip_margin, 0)
            ),
        )
    )
    if mask:
        alignment_cell = mask_alignment_cross(
            bottom_layer=trench_layer, top_layer=feeding_channel_layer
        )
    else:
        alignment_cell = alignment_cross(layer=trench_layer)
    if alignment_mark_position is None:
        alignment_mark_position = chip_area_corner[1] * 7 / 6
    if mask:
        vertical_alignment_spacing = np.array([0, alignment_mark_position * 2])
        main_cell.add(
            Reference(
                alignment_cell,
                -vertical_alignment_spacing / 2,
                columns=1,
                rows=2,
                spacing=vertical_alignment_spacing,
            )
        )
        horizontal_alignment_spacing = np.array([alignment_mark_position * 2, 0])
        main_cell.add(
            Reference(
                alignment_cell,
                -horizontal_alignment_spacing / 2,
                columns=2,
                rows=1,
                spacing=horizontal_alignment_spacing,
            )
        )
    else:
        alignment_spacing = np.array([0, alignment_mark_position * 2])
        main_cell.add(
            Reference(
                alignment_cell,
                -alignment_spacing / 2,
                columns=1,
                rows=2,
                spacing=alignment_spacing,
            )
        )
    if label:
        if mask:
            mask_label_padding = 30e2
            fc_label_position = corner - mask_label_padding
            trench_label_position = fc_label_position * np.array([-1, 1])
            main_cell.add(
                *text(
                    "FC layer\n" + label_right,
                    label_font_size,
                    position=fc_label_position,
                    horizontal_alignment="right",
                    layer=feeding_channel_layer,
                )
            )
            main_cell.add(
                *text(
                    "trench layer\n" + label_right,
                    label_font_size,
                    position=trench_label_position,
                    horizontal_alignment="left",
                    layer=trench_layer,
                )
            )
        else:
            label_right_x = (
                chip_area_corner[0]
                + np.abs(profilometry_bbox[1, 0] - profilometry_bbox[0, 0])
                + 2 * chip_margin
            )
            label_right_position = (label_right_x, 0)
            if label_right:
                main_cell.add(
                    *text(
                        label_right,
                        right_font_size,
                        position=label_right_position,
                        angle=np.pi / 2,
                        horizontal_alignment="left",
                        prerotate_alignment="center",
                        layer=feeding_channel_layer,
                    )
                )
            label_left_position = (-label_right_x, 0)
            if label_left:
                main_cell.add(
                    *text(
                        label_left,
                        left_font_size,
                        position=label_left_position,
                        angle=-np.pi / 2,
                        horizontal_alignment="right",
                        prerotate_alignment="center",
                        layer=feeding_channel_layer,
                    )
                )
    return main_cell


def chip(
    name,
    design_func=snake,
    ul_corner_label=None,
    ur_corner_label=None,
    ll_corner_label=None,
    lr_corner_label=None,
    corner_label_font_size=1500,
    corner_label_margin=400,
    feeding_channel_layer=FEEDING_CHANNEL_LAYER,
    trench_layer=TRENCH_LAYER,
    metadata=None,
    **kwargs,
):
    """Generates a GDS cell containing a chip. Keyword arguments are passed
    through to `design_func`, which is expected to return a GDS cell containing
    a chip design. This function adds an outline and a text label to that
    design.

    Parameters
    ----------
    name : str
        Human-readable label, will be written on the top edge of the chip.
    design_func : function, optional
        Function which will return a GDS cell containing the chip design.
    ul_corner_label : str, optional
        Label to put in the upper left corner of the chip.
    ur_corner_label : str, optional
        Label to put in the upper right corner of the chip.
    ll_corner_label : str, optional
        Label to put in the lower left corner of the chip.
    lr_corner_label : str, optional
        Label to put in the lower right corner of the chip.
    corner_label_font_size : float, optional
        Font size for corner labels (in µm).
    feeding_channel_layer : int, optional
        Feeding channel layer number.
    trench_layer : int, optional
        Trench layer number.
    metadata : Union[dict, None]
        If not None, add metadata to the dict.

    Returns
    -------
    chip_cell
        GDS cell containing the outline, text label, and chip design.
    """
    if "dims" not in kwargs:
        kwargs["dims"] = DEFAULT_DIMS
    kwargs["feeding_channel_layer"] = feeding_channel_layer
    kwargs["trench_layer"] = trench_layer
    dims = kwargs["dims"]
    design_cell, md = design_func(**{"name": name, **kwargs})
    if metadata is not None:
        metadata[name] = md  # TODO: this won't work with memoization!!!
    chip_cell = Cell(f"Chip-{name}")
    chip_cell.add(*outline(dims, layer=feeding_channel_layer))
    corner_x = dims[0] / 2 - corner_label_margin
    corner_y = dims[1] / 2 - corner_label_margin
    if ul_corner_label:
        chip_cell.add(
            *text(
                ul_corner_label,
                corner_label_font_size,
                (-corner_x, corner_y),
                horizontal_alignment="left",
                vertical_alignment="top",
                layer=feeding_channel_layer,
            )
        )
    if ur_corner_label:
        chip_cell.add(
            *text(
                ur_corner_label,
                corner_label_font_size,
                (corner_x, corner_y),
                horizontal_alignment="right",
                vertical_alignment="top",
                layer=feeding_channel_layer,
            )
        )
    if ll_corner_label:
        chip_cell.add(
            *text(
                ll_corner_label,
                corner_label_font_size,
                (-corner_x, -corner_y),
                horizontal_alignment="left",
                vertical_alignment="bottom",
                layer=feeding_channel_layer,
            )
        )
    if lr_corner_label:
        chip_cell.add(
            *text(
                lr_corner_label,
                corner_label_font_size,
                (corner_x, -corner_y),
                horizontal_alignment="right",
                vertical_alignment="bottom",
                layer=feeding_channel_layer,
            )
        )
    chip_cell.add(Reference(design_cell, (0, 0)))
    return chip_cell


def outline(dims, thickness=0.15e3, layer=FEEDING_CHANNEL_LAYER):
    """Generates a rectangular outline marking the edges of a chip.

    Parameters
    ----------
    dims : Tuple[float, float]
        A tuple of `(height, width)` specifying the dimensions of the chip footprint (in
        microns).
    thickness : float, optional
        Thickness (in microns) of the outline.
    layer : int, optional
        Layer number.

    Returns
    -------
    outline
        GDS cell containing the outline.
    """
    outline_inner = rectangle(-dims / 2, dims / 2)
    outline_outer = rectangle(-(dims + thickness) / 2, (dims + thickness) / 2)
    outline = boolean(outline_outer, outline_inner, "not", layer=layer)
    return outline
