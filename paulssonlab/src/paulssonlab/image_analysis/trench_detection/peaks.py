import warnings

import holoviews as hv
import numpy as np
import pandas as pd
import scipy.ndimage
import scipy.signal


def find_peaks(
    profile,
    prominence=0.2,
    min_distance=10,
    prominence_wlen=None,
    detrend_window=200,
    diagnostics=None,
):
    if prominence_wlen is None:
        prominence_wlen = 5 * min_distance
    minimum = scipy.ndimage.minimum_filter1d(profile, detrend_window, mode="constant")
    maximum = scipy.ndimage.maximum_filter1d(profile, detrend_window, mode="constant")
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            "invalid value encountered in true_divide",
        )
        profile_detrended = (profile - minimum) / (maximum - minimum)
    # TODO: which properties do we want to calculate?
    idxs, properties = scipy.signal.find_peaks(
        profile_detrended,
        distance=min_distance,
        height=(None, None),
        prominence=prominence,
        wlen=prominence_wlen,
        width=(None, None),
    )
    if diagnostics is not None:
        points = hv.Scatter((idxs, profile[idxs])).options(size=5, color="cyan")
        points_detrended = hv.Scatter((idxs, profile_detrended[idxs])).options(
            size=5, color="cyan"
        )
        diagnostics["profile_detrended"] = (
            hv.Curve(profile)
            * hv.Curve(minimum).opts(color="red")
            * hv.Curve(maximum).opts(color="green")
            * points
        )
        diagnostics["peaks"] = hv.Curve(profile_detrended) * points_detrended
    trench_info = pd.DataFrame(properties)
    return idxs, None, trench_info


def find_periodic_peaks(
    profile,
    pitch=None,
    refine=5,
    nfft=2**14,
    smooth_offset=4,
    num_offset_points=200,
    min_period=50,
    diagnostics=None,
):
    if pitch is None or diagnostics is not None:
        freqs, spectrum = scipy.signal.periodogram(
            profile, window="hann", nfft=nfft, scaling="spectrum"
        )
        if min_period:
            spectrum[:min_period] = 0
    if pitch is None:
        pitch_idx = spectrum.argmax()
        pitch = 1 / freqs[pitch_idx]
    if diagnostics is not None:
        diagnostics["pitch"] = pitch
        spectrum_plot = hv.Curve((1 / freqs, spectrum))
        spectrum_plot *= hv.VLine(pitch).options(color="red")
        diagnostics["spectrum"] = spectrum_plot
    offsets = np.linspace(0, pitch, num_offset_points, endpoint=False)
    # always leave a period of length `pitch` so we can add offsets
    # even when we could've fit another period
    # TODO: this sometimes results in one trench undetected!
    offset_idxs = (
        np.arange(0, len(profile) - pitch, pitch) + offsets[:, np.newaxis]
    ).astype(np.int_)
    offset_objective = profile[offset_idxs].sum(axis=1)
    if smooth_offset:
        offset_objective_smoothed = scipy.ndimage.filters.gaussian_filter1d(
            offset_objective, smooth_offset
        )
    else:
        offset_objective_smoothed = offset_objective
    offset_idx = offset_objective_smoothed.argmax()
    offset = offsets[offset_idx]
    idxs = offset_idxs[offset_idx]
    if diagnostics is not None:
        diagnostics["offset"] = offset
        offset_plot = hv.Curve((offsets, offset_objective))
        if smooth_offset:
            offset_plot *= hv.Curve((offsets, offset_objective_smoothed)).options(
                color="cyan"
            )
        offset_plot *= hv.VLine(offset).options(color="red")
        diagnostics["offsets"] = offset_plot
    if not refine:
        refined_idxs = idxs
        trench_info = {}
    else:
        idx_start = np.clip(idxs - refine, 0, len(profile))
        idx_end = np.clip(idxs + refine, 0, len(profile))
        shifts = np.array(
            [
                profile[idx_start[i] : idx_end[i]].argmax() - (idxs[i] - idx_start[i])
                for i in range(len(idxs))
            ]
        )
        shifts = np.where(profile[idxs] != profile[idxs + shifts], shifts, 0)
        refined_idxs = idxs + shifts
        trench_info = pd.DataFrame({"shift": shifts})
        if diagnostics is not None:
            periodic_points = hv.Scatter((idxs, profile[idxs])).options(
                size=5, color="red"
            )
            refined_points = hv.Scatter((refined_idxs, profile[refined_idxs])).options(
                size=3, color="cyan"
            )
            diagnostics["refined_points"] = (
                hv.Curve(profile) * periodic_points * refined_points
            )
    prominence_data = scipy.signal.peak_prominences(profile, refined_idxs)
    width_data = scipy.signal.peak_widths(
        profile, refined_idxs, prominence_data=prominence_data
    )
    trench_info = {
        **trench_info,
        **dict(zip(("prominences", "left_bases", "right_bases"), prominence_data)),
        **dict(zip(("widths", "width_heights", "left_ips", "right_ips"), width_data)),
    }
    trench_info = pd.DataFrame(trench_info)
    info = dict(pitch=pitch, offset=offset)
    return refined_idxs, info, trench_info
