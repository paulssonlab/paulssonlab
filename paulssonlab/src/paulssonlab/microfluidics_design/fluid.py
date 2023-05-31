from collections.abc import Mapping

import numpy as np
import pint

from paulssonlab.microfluidics_design.util import strip_units


def _ensure_array(val, size):
    if np.isscalar(val):
        return np.full(size, val)
    else:
        return val


def entrance_length(height, width, flow_rate, kinematic_viscosity):
    """Calculates entrance length (distance to fully-developed laminar flow)
    using equation 8 from [^1].

    [^1]: Ahmad, T., & Hassan, I. (2010). Experimental analysis of microchannel entrance length characteristics using microparticle image velocimetry. Journal of Fluids Engineering, 132(4).

    Parameters
    ----------
    height : float
        Height of channel.
    width : float
        Width of channel.
    flow_rate : float
        Volumetric flow rate.
    kinematic_viscosity : float
        Kinematic viscosity of the fluid.

    Returns
    -------
    L_e
        Entrance length.
    """
    h = height
    w = width
    Q = flow_rate
    nu = kinematic_viscosity
    # hydraulic diameter
    D_h = 2 * w * h / (w + h)
    # Reynolds number
    Re = Q * D_h / ((w * h) * nu)
    # entry length
    L_e = D_h * (0.55 / (0.13 * Re + 1) + 0.065 * Re)
    return L_e


def resistance(height, width):
    """Calculates hydraulic resistance (per length per dynamic viscosity) along
    a rectangular channel. Multiply by channel length and dynamic viscosity to
    get hydraulic resistance with the correct units (i.e., obeys `Q=delta_P/R`
    where `Q` is the volumetric flow rate and `delta_P` is the pressure drop
    along the channel).

    To derive, start with the mean velocity
    `v = delta_P * (h/2)**2 * (1/3 - 64*eps/np.pi**5*np.tanh(np.pi/2/eps))`
    (from [^1]) where `eps = h/w`, `h <= w`, and `c=h/2`. Then use that
    `Q = v * w * h` and `R = delta_P / Q`. For the simplified expression, see [^2].

    [^1]: Bahrami, M., Yovanovich, M. M., & Culham, J. R. (2006). Pressure drop of fully-developed, laminar flow in microchannels of arbitrary cross-section.
    [^2]: https://www.elveflow.com/microfluidic-reviews/microfluidic-flow-control/flow-control-in-microfluidics-device/

    Parameters
    ----------
    height : float
        Height of channel.
    width : float
        Width of channel.

    Returns
    -------
    resistance
    """
    h = np.minimum(height, width)
    w = np.maximum(height, width)
    Q = h**3 * w / 12 * (1 - 0.63 * h / w * np.tanh(1.57 * w / h))
    return 1 / Q


def bend_resistance(height, width, inner_radius, angle):
    return 0


def ladder_flow_rates(R_snake, R_left, R_right, N=None):
    for ary, inc in ((R_snake, 0), (R_left, 1), (R_right, 1)):
        if not np.isscalar(ary):
            if N is None:
                N = len(ary) + inc
            elif N != len(ary) + inc:
                raise ValueError("got conflicting ladder sizes")
    if N is None:
        raise ValueError("ladder size must be specified if resistances are scalars")
    R_snake = _ensure_array(R_snake, N)
    R_left = _ensure_array(R_left, N - 1)
    R_right = _ensure_array(R_right, N - 1)
    A = np.zeros((N, N))
    b = np.zeros(N)
    b[0] = 1
    A[0, :] = 1
    for i in range(1, N):
        A[i, : i - 1] = R_left[i - 1]
        A[i, i - 1] = -(R_right[i - 1] + R_snake[i - 1])
        A[i, i] = R_left[i - 1] + R_snake[i]
        A[i, i + 1 :] = -R_right[i - 1]
    q = np.linalg.solve(A, b)
    return q


def manifold_flow_rates(
    feeding_channel_height=None,
    snake_split=None,
    manifold_split=None,
    input_info=None,
    **kwargs,
):
    snake_split_cum = np.concatenate(((0,), np.cumsum(snake_split)))
    manifold_split_cum = np.concatenate(((0,), np.cumsum(manifold_split)))
    results = {}
    for input_name, input_md in input_info.items():
        lane_length = input_md["lane_length"]
        manifold_width = input_md["manifold_width"]
        feeding_channel_width = input_md["feeding_channel_width"]
        inner_snake_bend_radius = input_md["inner_snake_bend_radius"]
        left_port_lane_ys = input_md["left_port_lane_ys"]
        right_port_lane_ys = input_md["right_port_lane_ys"]
        lane_ys = input_md["lane_ys"]
        lanes_per_snake = input_md["lanes_per_snake"]
        left_segment_lengths = -np.diff(left_port_lane_ys)
        right_segment_lengths = -np.diff(right_port_lane_ys)
        R_snake = lanes_per_snake * lane_length * resistance(
            feeding_channel_height, feeding_channel_width
        ) + (lanes_per_snake - 1) * bend_resistance(
            feeding_channel_height,
            feeding_channel_width,
            inner_snake_bend_radius,
            np.pi,
        )
        R_left = left_segment_lengths * resistance(
            feeding_channel_height, manifold_width
        )
        R_right = right_segment_lengths * resistance(
            feeding_channel_height, manifold_width
        )
        q = ladder_flow_rates(
            strip_units(R_snake), strip_units(R_left), strip_units(R_right)
        )
        # TODO: is this correct?
        R = ((1 - np.cumsum(q))[:-1] * R_left).sum() + q[-1] * R_snake[-1]
        results[input_name] = {
            "R_snake": R_snake,
            "R_left": R_left,
            "R_right": R_right,
            "q": q,
            "R": R,
        }
    return results


def _metadata_to_units(metadata, length_unit):
    metadata = {
        **metadata,
        **{
            k: metadata[k] * length_unit
            for k in (
                "lane_length",
                "feeding_channel_width",
                "manifold_width",
                "inner_snake_bend_radius",
                "lane_ys",
                "left_port_lane_ys",
                "right_port_lane_ys",
            )
            if k in metadata and metadata[k] is not None
        },
        **{
            k: _metadata_to_units(metadata[k], length_unit)
            for k in metadata
            if isinstance(metadata[k], Mapping)
        },
    }
    return metadata


def snake_flow(
    feeding_channel_height,
    metadata,
    length_unit=None,
):
    """Summary.

    Parameters
    ----------
    feeding_channel_height : TYPE
        Description
    metadata : TYPE
        Description
    length_unit : None, optional
        Description

    Returns
    -------
    TYPE
        Description
    """
    if length_unit is None:
        length_unit = 1
    metadata = _metadata_to_units(metadata, length_unit)
    res = {}
    if "manifold_split" in metadata:
        flow_info = manifold_flow_rates(feeding_channel_height, **metadata)
        for input_name, input_flow in flow_info.items():
            q = input_flow["q"]
            flow_nonuniformity = q.min() / q.max()
            res[input_name] = {**input_flow, "flow_nonuniformity": flow_nonuniformity}

    else:
        for input_name, input_md in metadata["input_info"].items():
            R = input_md["snake_length"] * resistance(
                feeding_channel_height, input_md["feeding_channel_width"]
            )
            res[input_name] = {"R": R}
    return res
