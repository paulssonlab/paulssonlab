#!/usr/bin/env python
import numpy as np
import gdspy as g
from gdspy import CellReference, CellArray, Rectangle
from geometry import (
    MAX_POINTS,
    ROUND_POINTS,
    Cell,
    Round,
    cross,
    qr_target,
    boolean,
    mirror,
    mirror_refs,
    align,
    align_refs,
    flatten_or_merge,
)
from text import Text as _Text
from util import make_odd, memoize
import hamming
import bitarray.util
import click
from functools import partial
from cytoolz import compose
from itertools import product
import numbers
import shortuuid

REFERENCE_LAYER = 0
TRENCH_LAYER = 1
FEEDING_CHANNEL_LAYER = 2
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
    if prerotate_alignment is not None:
        objs = align_refs(objs, position=(0, 0), alignment=prerotate_alignment)
    if angle == 0:
        pass
    elif angle == -np.pi / 2:
        for ref in objs:
            ref.rotation += 90
            ref.origin[:] = ref.origin[::-1] * np.array([-1, 1])
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

# TODO
# fix duplicate chip memoization

# @memoize
def sampler_snake(trench_width=[1.5], trench_length=[35], **kwargs):
    """Generates a snake with different trench widths and lengths in a single
    chip.

    Parameters
    ----------
    trench_width : List[float], optional
        Trench widths (in microns).
    trench_length : List[float], optional
        Trench lengths (in microns).
    **kwargs
        Other keyword arguments are passed through to `snake`.
    """
    raise NotImplementedError


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


# @memoize # TODO!!!!
def manifold_snake(
    dims=DEFAULT_DIMS,
    split=None,
    add_remainder=False,
    lanes_per_snake=1,
    num_manifolds=1,
    manifold_width=200,
    manifold_input_margin=2e3,
    manifold_bend_margin=0.2e3,
    manifold_margin=200,
    manifold_bend_radius=200,
    manifold_round_radius=True,
    manifold_input_style="bend-out",
    border_margin=600,
    top_margin=1.5e3,
    bottom_margin=1e3,
    port_margin=0.5e3,
    trench_width=1.5,
    trench_length=35,
    trench_fc_overlap=None,
    trench_margin=0.5e3,
    trench_gap=20,
    trench_spacing=2,
    feeding_channel_width=90,
    port_radius=200,
    registration_marks=False,
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
    tick_text_size=None,
    tick_labels=False,
    draw_trenches=True,
    flatten_feeding_channel=False,
    merge_feeding_channel=True,
    feeding_channel_layer=FEEDING_CHANNEL_LAYER,
    trench_layer=TRENCH_LAYER,
    label=None,
):
    if manifold_input_style not in ("u-turn", "bend-out", "bend-in"):
        raise NotImplementedError
    if label is None:
        label = get_uuid()
    if trench_fc_overlap is None:
        trench_fc_overlap = min(trench_length, feeding_channel_width / 3)
    if tick_text_size is None:
        tick_text_size = tick_length * 2
    if manifold_round_radius is True:
        manifold_round_radius = feeding_channel_width
    # define some dimensions
    effective_trench_length = trench_length + trench_gap / 2
    inner_snake_bend_radius = effective_trench_length
    outer_snake_bend_radius = feeding_channel_width + inner_snake_bend_radius
    if manifold_input_style == "u-turn":
        horizontal_margin = (
            port_margin
            + port_radius
            + manifold_width / 2
            + 2 * manifold_bend_radius
            + manifold_width
            + manifold_margin
        )
    elif manifold_input_style == "bend-out":
        horizontal_margin = (
            port_margin
            + 2 * port_radius
            + manifold_input_margin
            + manifold_bend_radius
            + manifold_width
            + manifold_margin
        )
    elif manifold_input_style == "bend-in":
        horizontal_margin = border_margin + manifold_width + manifold_margin
        # TODO: rename manifold_margin to something more clear
        # TODO: combine border/port margin?
        # TODO: this should be +=, so that top/bottom_margin is extra
        # TODO: top_margin should incorporate label size
        top_margin = (
            port_margin
            + port_radius
            + manifold_width / 2
            + manifold_bend_radius
            + manifold_bend_margin
        )
        bottom_margin = top_margin
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
        if split:
            raise ValueError("cannot specify both split and lanes_per_snake")
        split = (lanes_per_snake,) * int(max_lanes // lanes_per_snake)
        if add_remainder:
            remainder = max_lanes % lanes_per_snake
            if remainder % 2 == 0:
                remainder -= 1
            if remainder > 0:
                split += (remainder,)
    split = _compute_lane_split(split, max_lanes)
    num_lanes = np.sum(split)
    # manifold split
    split_cum = np.concatenate(((0,), np.cumsum(split)))
    left_port_lanes = split_cum[1:] - 1
    right_port_lanes = split_cum[:-1]
    lanes_per_input = int(len(split) // num_manifolds)
    manifold_split = (lanes_per_input,) * num_manifolds
    remainder = len(split) % num_manifolds
    if remainder > 0:
        manifold_split = (*manifold_split[:-1], manifold_split[-1] + remainder)
    manifold_split_cum = np.concatenate(((0,), np.cumsum(manifold_split)))
    # define trench parameters
    trench_xs = np.arange(
        -(lane_fc_dims[0] - trench_margin - trench_width) / 2,
        (lane_fc_dims[0] - trench_margin - trench_width) / 2,
        trench_width + trench_spacing,
    )
    trenches_per_set = len(trench_xs)
    num_trenches = trenches_per_set * 2 * num_lanes
    # feeding channel
    snake_cell = Cell(f"Snake-{label}")
    snake_fc_cell, lane_ys = _snake_feeding_channel(
        lane_fc_dims,
        effective_trench_length,
        outer_snake_bend_radius,
        0,
        split,
        label,
        layer=feeding_channel_layer,
    )
    # manifolds
    if manifold_round_radius:
        rounded_corner = Cell(f"Snake-round-{label}")
        rounded_curve = g.Curve(0, 0)
        rounded_curve.v(manifold_round_radius)
        rounded_curve.arc(manifold_round_radius, 0, -1 / 2 * np.pi, 0)
        rounded_corner.add(
            g.Polygon(rounded_curve.get_points(), layer=feeding_channel_layer)
        )
    for idx in range(len(manifold_split_cum) - 1):
        for flip in (1, -1):
            snake_manifold_cell = Cell(
                f"Snake-{'left' if flip == 1 else 'right'}-{idx}-{label}"
            )
            port_lanes = right_port_lanes if flip == -1 else left_port_lanes
            manifold_lane_ys = lane_ys[
                port_lanes[manifold_split_cum[idx] : manifold_split_cum[idx + 1]]
            ][::-flip]
            manifold_input_bend_y = manifold_lane_ys[0] - flip * (
                feeding_channel_width / 2 + manifold_bend_margin
            )
            if manifold_input_style == "u-turn":
                port_x = -dims[0] / 2 + port_margin + port_radius
                manifold_input_bend_x = (
                    port_x + manifold_width / 2 + manifold_bend_radius
                )
                manifold_left_x = manifold_input_bend_x + manifold_bend_radius
                port_y = manifold_input_bend_y + flip * (
                    manifold_input_margin + port_radius
                )
                manifold_bend_angles = (0, -flip * np.pi)
            elif manifold_input_style == "bend-out":
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
            snake_manifold_cell.add(
                Round(
                    (-flip * port_x, port_y), port_radius, layer=feeding_channel_layer
                )
            )
            snake_manifold_cell.add(
                Round(
                    (-flip * manifold_input_bend_x, manifold_input_bend_y),
                    manifold_bend_radius + manifold_width,
                    inner_radius=manifold_bend_radius,
                    initial_angle=manifold_bend_angles[0],
                    final_angle=manifold_bend_angles[1],
                    layer=feeding_channel_layer,
                )
            )
            if manifold_input_style in ("bend-out", "bend-in"):
                snake_manifold_cell.add(
                    Rectangle(
                        (-flip * port_x, port_y + manifold_width / 2),
                        (-flip * manifold_input_bend_x, port_y - manifold_width / 2),
                        layer=feeding_channel_layer,
                    )
                )
            elif manifold_input_style == "u-turn":
                snake_manifold_cell.add(
                    Rectangle(
                        (-flip * (port_x - manifold_width / 2), manifold_input_bend_y),
                        (-flip * (port_x + manifold_width / 2), port_y),
                        layer=feeding_channel_layer,
                    )
                )
            # manifold bend margin
            snake_manifold_cell.add(
                Rectangle(
                    (-flip * manifold_left_x, manifold_input_bend_y),
                    (-flip * (manifold_left_x + manifold_width), manifold_top_y),
                    layer=feeding_channel_layer,
                )
            )
            # manifold
            snake_manifold_cell.add(
                Rectangle(
                    (-flip * manifold_left_x, manifold_top_y),
                    (-flip * (manifold_left_x + manifold_width), manifold_taper_y),
                    layer=feeding_channel_layer,
                )
            )
            if manifold_round_radius:
                for y in manifold_lane_ys[:-1]:
                    snake_manifold_cell.add(
                        CellReference(
                            rounded_corner,
                            (
                                -flip * (manifold_left_x + manifold_width),
                                y + flip * feeding_channel_width / 2,
                            ),
                            rotation=90 * (3 + flip),
                        )
                    )
                for y in manifold_lane_ys:
                    snake_manifold_cell.add(
                        CellReference(
                            rounded_corner,
                            (
                                -flip * (manifold_left_x + manifold_width),
                                y - flip * feeding_channel_width / 2,
                            ),
                            rotation=90 * (4 + flip),
                        )
                    )
            for y in manifold_lane_ys:
                snake_manifold_cell.add(
                    g.Rectangle(
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
            snake_fc_cell.add(
                Round(
                    (-flip * (manifold_left_x + manifold_width), manifold_taper_y),
                    manifold_width,
                    initial_angle=manifold_taper_angle,
                    final_angle=manifold_taper_angle + flip * np.pi,
                    layer=feeding_channel_layer,
                )
            )
            snake_fc_cell.add(CellReference(snake_manifold_cell, (0, 0)))
    flatten_or_merge(
        snake_fc_cell,
        flatten=flatten_feeding_channel,
        merge=merge_feeding_channel,
        layer=feeding_channel_layer,
    )
    y_offset = (bottom_margin - top_margin) / 2
    snake_cell.add(CellReference(snake_fc_cell, (0, y_offset)))
    # trenches
    if draw_trenches:
        snake_trenches_cell = _snake_trenches(
            trench_width=trench_width,
            trench_spacing=trench_spacing,
            trench_length=trench_length,
            trench_gap=trench_gap,
            trench_fc_overlap=trench_fc_overlap,
            feeding_channel_width=feeding_channel_width,
            registration_marks=registration_marks,
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
            tick_text_size=tick_text_size,
            trench_xs=trench_xs,
            lane_ys=lane_ys,
            label=label,
            layer=TRENCH_LAYER,
        )
        snake_cell.add(CellReference(snake_trenches_cell, (0, y_offset)))
    metadata = {
        k: v
        for k, v in locals().items()
        if k
        in (
            "num_lanes",
            "trenches_per_set",
            "num_trenches",
            "feeding_channel_width",
            "trench_gap",
            "trench_length",
            "manifold_width",
            "outer_snake_bend_radius",
            "lane_ys",
            "split",
            "split_cum",
            "manifold_split",
            "manifold_split_cum",
            "left_port_lanes",
            "right_port_lanes",
        )
    }
    metadata["lane_length"] = lane_fc_dims[0]
    metadata["lane_with_trenches_length"] = trench_xs[-1] - trench_xs[0] + trench_width
    return snake_cell, metadata


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
    trench_gap=20,
    gap_lanes=0,
    trench_spacing=2,
    feeding_channel_width=90,
    port_radius=200,
    registration_marks=False,
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
    tick_text_size=None,
    tick_labels=False,
    draw_trenches=True,
    flatten_feeding_channel=False,
    merge_feeding_channel=True,
    feeding_channel_layer=FEEDING_CHANNEL_LAYER,
    trench_layer=TRENCH_LAYER,
    label=None,
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
    horizontal_margin : float, optional
        Empty space (in microns) between the outside bends of the snake and the chip
        border.
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
    tick_text_size : None, optional
        Font size (in microns) of the trench tick numbering.
    tick_labels : bool, optional
        Whether to label every tick mark with the trench number.
    draw_trenches : bool, optional
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
    label : str, optional
        Unique identifier used for cell names. If None is given, an alphanumeric UUID
        will be generated.

    Returned
    ------------------
    snake_cell
        GDS cell containing the chip design.
    metadata : dict
        Dict containing the design parameters.
    """
    if label is None:
        label = get_uuid()
    if port_margin is None:
        port_margin = horizontal_margin / 2
    if trench_fc_overlap is None:
        trench_fc_overlap = min(trench_length, feeding_channel_width / 3)
    if tick_text_size is None:
        tick_text_size = tick_length * 2
    effective_trench_length = trench_length + trench_gap / 2
    inner_snake_bend_radius = effective_trench_length
    outer_snake_bend_radius = feeding_channel_width + inner_snake_bend_radius
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
            "trench_gap",
            "trench_length",
        )
    }
    metadata["lane_length"] = lane_fc_dims[0]
    metadata["lane_with_trenches_length"] = trench_xs[-1] - trench_xs[0] + trench_width
    metadata["max_snake_length"] = max(split) * lane_fc_dims[0]
    snake_cell = Cell(f"Snake-{label}")
    port_offset = dims[0] / 2 - lane_fc_dims[0] / 2 - port_radius - port_margin
    snake_fc_cell, lane_ys = _snake_feeding_channel(
        lane_fc_dims,
        effective_trench_length,
        port_offset,
        port_radius,
        split,
        label,
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
    snake_cell.add(CellReference(snake_fc_cell, (0, y_offset)))
    if draw_trenches:
        snake_trenches_cell = _snake_trenches(
            trench_width=trench_width,
            trench_spacing=trench_spacing,
            trench_length=trench_length,
            trench_gap=trench_gap,
            trench_fc_overlap=trench_fc_overlap,
            feeding_channel_width=feeding_channel_width,
            registration_marks=registration_marks,
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
            tick_text_size=tick_text_size,
            trench_xs=trench_xs,
            lane_ys=lane_ys,
            label=label,
            layer=TRENCH_LAYER,
        )
        snake_cell.add(CellReference(snake_trenches_cell, (0, y_offset)))
    return snake_cell, metadata


def _snake_feeding_channel(
    lane_fc_dims,
    effective_trench_length,
    port_offset,
    port_radius,
    split,
    label,
    gap_lanes=0,
    layer=FEEDING_CHANNEL_LAYER,
    flatten_feeding_channel=False,
    merge_feeding_channel=True,
):
    feeding_channel_width = lane_fc_dims[1]
    lane_height = feeding_channel_width + 2 * effective_trench_length
    inner_snake_bend_radius = effective_trench_length
    outer_snake_bend_radius = feeding_channel_width + inner_snake_bend_radius
    lane_cell = Cell(f"Lane-{label}")
    lane_fc = Rectangle(-lane_fc_dims / 2, lane_fc_dims / 2, layer=layer)
    lane_cell.add(lane_fc)
    bend_cell = Cell(f"Feeding Channel Bend-{label}")
    bend = Round(
        (0, 0),
        outer_snake_bend_radius,
        inner_radius=inner_snake_bend_radius,
        initial_angle=3 / 2 * np.pi,
        final_angle=5 / 2 * np.pi,
        layer=layer,
    )
    bend_cell.add(bend)
    port_cell = Cell(f"Feeding Channel Port-{label}")
    if port_radius:
        port = Round((port_offset, 0), port_radius, layer=layer)
        port_cell.add(port)
    port_fc = Rectangle(
        (0, -feeding_channel_width / 2),
        (port_offset, feeding_channel_width / 2),
        layer=layer,
    )
    port_cell.add(port_fc)
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
    snake_fc_cell = Cell(f"Snake Feeding Channel-{label}")
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
    return snake_fc_cell, lane_ys


@memoize  # TODO!!!
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
        zero_symbol = g.Polygon([(1 / 2, -1 / 2), (1 / 2, 1 / 2), (-1 / 2, -1 / 2)])
    elif zero_symbol is not False:
        zero_symbol = gdspy.copy(zero_symbol)
    if one_symbol is None:
        one_symbol = g.Rectangle((-1 / 2, -1 / 2), (1 / 2, 1 / 2))
    elif one_symbol is not False:
        one_symbol = gdspy.copy(one_symbol)
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
                    cell.add(g.copy(one_symbol, x, y))
            else:
                if zero_symbol is not False:
                    cell.add(g.copy(zero_symbol, x, y))
    return cell


def _snake_trenches(
    trench_width,
    trench_spacing,
    trench_length,
    trench_gap,
    trench_fc_overlap,
    feeding_channel_width,
    registration_marks,
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
    tick_text_size,
    trench_xs,
    lane_ys,
    label,
    layer=TRENCH_LAYER,
):
    if ticks and registration_marks:
        raise ValueError("cannot draw both ticks and registration marks")
    lane_gap_offset_y = feeding_channel_width / 2 + trench_length + trench_gap / 2
    mark_pitch = mark_size + mark_spacing
    column_barcode_margin = (
        4 * mark_size + 3 * mark_spacing
    ) / 2 + mark_spacing  # TODO: couple to qr_target arguments, below
    row_barcode_margin = column_barcode_margin + barcode_columns * mark_pitch
    chip_barcode_margin = row_barcode_margin + barcode_columns * mark_pitch
    trenches_per_set = len(trench_xs)
    snake_trenches_cell = Cell(f"Snake Trenches-{label}")
    trench_cell = Cell(f"Trench-{label}")
    trench_cell.add(
        Rectangle(
            (-trench_width / 2, -trench_fc_overlap),
            (trench_width / 2, trench_length),
            layer=layer,
        )
    )
    tick_cell = Cell(f"Tick-{label}")
    if ticks:
        tick_cell.add(
            Rectangle(
                (-trench_width / 2, tick_margin - trench_gap / 2),
                (trench_width / 2, tick_margin - trench_gap / 2 + tick_length),
                layer=layer,
            )
        )
    elif registration_marks:
        tick_cell.add(
            qr_target(
                mark_size, mark_spacing, 2 * mark_size + mark_spacing, layer=layer
            )
        )
    tick_xs = trench_xs[::tick_period]
    num_ticks = len(tick_xs)
    if registration_marks:
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
                    CellArray(
                        chip_barcode,
                        num_ticks,
                        len(lane_ys),
                        (
                            tick_period * (trench_width + trench_spacing),
                            lane_ys_diff[0],
                        ),
                        (
                            trench_xs[0] + chip_barcode_margin,
                            lane_ys[0] + lane_gap_offset_y,
                        ),
                    )
                )
            else:
                for y in lane_ys:
                    snake_trenches_cell.add(
                        CellArray(
                            chip_barcode,
                            num_ticks,
                            1,
                            (tick_period * (trench_width + trench_spacing), 0),
                            (trench_xs[0] + row_barcode_margin, y + lane_gap_offset_y),
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
                    CellArray(
                        column_barcode,
                        1,
                        len(lane_ys),
                        (0, lane_ys_diff[0]),
                        (x + column_barcode_margin, lane_ys[0] + lane_gap_offset_y),
                    )
                )
            else:
                for y in lane_ys:
                    snake_trenches_cell.add(
                        CellReference(
                            column_barcode,
                            (x + column_barcode_margin, y + lane_gap_offset_y),
                        )
                    )
    for lane_idx, y in enumerate(lane_ys):
        snake_trenches_cell.add(
            CellArray(
                trench_cell,
                trenches_per_set,
                1,
                (trench_width + trench_spacing, 0),
                (trench_xs[0], y + feeding_channel_width / 2),
            )
        )
        snake_trenches_cell.add(
            CellArray(
                trench_cell,
                trenches_per_set,
                1,
                (trench_width + trench_spacing, 0),
                (trench_xs[0], y - feeding_channel_width / 2),
                x_reflection=True,
            )
        )
        if ticks or registration_marks:
            snake_trenches_cell.add(
                CellArray(
                    tick_cell,
                    num_ticks,
                    1,
                    (tick_period * (trench_width + trench_spacing), 0),
                    (trench_xs[0], y + lane_gap_offset_y),
                )
            )
            if registration_marks:
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
                    CellArray(
                        row_barcode,
                        num_ticks,
                        1,
                        (tick_period * (trench_width + trench_spacing), 0),
                        (trench_xs[0] + row_barcode_margin, y + lane_gap_offset_y),
                    )
                )
        if tick_labels:
            for tick_idx, x in enumerate(tick_xs):
                tick_idx = (
                    lane_idx * 2 * trenches_per_set + 2 * tick_idx * tick_period + 1
                )
                snake_trenches_cell.add(
                    Text(
                        str(tick_idx),
                        tick_text_size,
                        (
                            x + 2 * trench_width,
                            y + feeding_channel_width / 2 + trench_length + tick_margin,
                        ),
                        layer=layer,
                    )
                )
    return snake_trenches_cell


# TODO: implement using bbox-aware auto-gridding helper
@memoize
def profilometry_marks(
    dims=np.array([150, 75]),
    columns=3,
    rows=3,
    layers=(FEEDING_CHANNEL_LAYER, TRENCH_LAYER),
    text=True,
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
    text : bool, optional
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
        Rectangle((-3 * length, -length), (-2 * length, length), layer=layer)
    )
    alignment_cell.add(
        Rectangle((2 * length, -length), (3 * length, length), layer=layer)
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
    outer_box = Rectangle(-(box_corner + thickness), box_corner + thickness)
    inner_box = Rectangle(-box_corner, box_corner)
    box = boolean(outer_box, inner_box, "not", layer=top_layer)
    cross_box_corner = np.array([length, length])
    base_cross_box = Rectangle(-cross_box_corner, cross_box_corner, layer=bottom_layer)
    base_cross = cross(length, thickness, layer=top_layer)
    base_cross_box = boolean(base_cross_box, base_cross, "not", layer=bottom_layer)
    alignment_cell.add(base_cross_box)
    alignment_cell.add(base_cross)
    for offset_unit_vector in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        offset_vector = np.array(offset_unit_vector) * offset_unit
        for num_offsets in range(1, num_crosses + 1):
            offset = offset_vector * num_offsets
            alignment_cell.add(g.copy(base_cross_box, *offset))
            alignment_cell.add(g.copy(base_cross, *offset))
    alignment_cell.add(box)
    return alignment_cell


def wafer(
    chips,
    text_right=None,
    text_left=None,
    diameter=76.2e3,
    chip_dims=None,
    chip_margin=0.5e3,
    alignment_mark_position=None,
    alignment_text_size=1000,
    right_text_size=2000,
    left_text_size=1300,
    text=True,
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
    text_right : str
        Text to write on right side of wafer.
    text_left : str
        Text to write on left side of wafer.
    diameter : float, optional
        Diameter of wafer (in microns).
    chip_dims : Tuple[float, float], optional
        Description
    chip_margin : float, optional
        Description
    alignment_mark_position : float, optional
        Distance (in microns) between both alignment marks.
    alignment_text_size : float, optional
        Size (in microns) of label text for alignment marks.
    right_text_size : float, optional
        Size (in microns) of text.
    left_text_size : float, optional
        Size (in microns) of text.
    text : bool, optional
        If true, label the wafer with `text_right` and `text_left`. If mask is true, additionally label the
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
    square = Rectangle(-square_corner, square_corner, layer=reference_layer)
    circle = Round((0, 0), diameter / 2)
    wafer_outline = boolean(square, circle, "not", layer=reference_layer)
    main_cell.add(wafer_outline)
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
        main_cell.add(CellReference(chip, (x, y)))
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
        layers=(feeding_channel_layer, trench_layer), text=text
    )
    profilometry_bbox = profilometry_cell.get_bounding_box()
    profilometry_spacing = np.array(
        [chip_area[0] + 2 * np.abs(profilometry_bbox[:, 0]).max(), 0]
    )
    main_cell.add(
        CellReference(
            profilometry_cell,
            np.array(
                (-chip_area[0] / 2 - profilometry_bbox[:, 0].max() - chip_margin, 0)
            ),
        )
    )
    main_cell.add(
        CellReference(
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
            CellArray(
                alignment_cell,
                1,
                2,
                vertical_alignment_spacing,
                -vertical_alignment_spacing / 2,
            )
        )
        horizontal_alignment_spacing = np.array([alignment_mark_position * 2, 0])
        main_cell.add(
            CellArray(
                alignment_cell,
                2,
                1,
                horizontal_alignment_spacing,
                -horizontal_alignment_spacing / 2,
            )
        )
    else:
        alignment_spacing = np.array([0, alignment_mark_position * 2])
        main_cell.add(
            CellArray(alignment_cell, 1, 2, alignment_spacing, -alignment_spacing / 2)
        )
    if text:
        if mask:
            mask_text_padding = 30e2
            fc_text_position = corner - mask_text_padding
            trench_text_position = fc_text_position * np.array([-1, 1])
            main_cell.add(
                Text(
                    "FC layer\n" + text_right,
                    label_text_size,
                    position=fc_text_position,
                    alignment="right",
                    layer=feeding_channel_layer,
                )
            )
            main_cell.add(
                Text(
                    "trench layer\n" + text_right,
                    label_text_size,
                    position=trench_text_position,
                    alignment="left",
                    layer=trench_layer,
                )
            )
        else:
            text_right_x = (
                chip_area_corner[0]
                + np.abs(profilometry_bbox[1, 0] - profilometry_bbox[0, 0])
                + 2 * chip_margin
            )
            text_right_position = (text_right_x, 0)
            if text_right:
                main_cell.add(
                    Text(
                        text_right,
                        right_text_size,
                        position=text_right_position,
                        angle=np.pi / 2,
                        alignment="left",
                        prerotate_alignment="centered",
                        layer=feeding_channel_layer,
                    )
                )
            text_left_position = (-text_right_x, 0)
            if text_left:
                main_cell.add(
                    Text(
                        text_left,
                        left_text_size,
                        position=text_left_position,
                        angle=-np.pi / 2,
                        alignment="right",
                        prerotate_alignment="centered",
                        layer=feeding_channel_layer,
                    )
                )
    return main_cell


# @memoize # TODO!!!
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
    label_text_size : float, optional
        Text size (in microns) of the label text.
    text : bool, optional
        If True, a label is written along the top edge of the chip.
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
    design_cell, md = design_func(**{"label": name, **kwargs})
    if metadata is not None:
        metadata[name] = md  # TODO: this won't work with memoization!!!
    chip_cell = Cell(f"Chip-{name}")
    chip_cell.add(outline(dims, layer=feeding_channel_layer))
    chip_cell.add(CellReference(design_cell, (0, 0)))
    text_position = (0, dims[1] / 2 - 1.1 * label_text_size)
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
    outline_inner = Rectangle(-dims / 2, dims / 2)
    outline_outer = Rectangle(-(dims + thickness) / 2, (dims + thickness) / 2)
    outline = boolean(outline_outer, outline_inner, "not", layer=layer)
    return outline


@click.group()
def cli():
    pass


if __name__ == "__main__":
    cli()
