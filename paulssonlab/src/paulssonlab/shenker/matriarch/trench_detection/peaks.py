import numpy as np
import pandas as pd
import holoviews as hv
import scipy.signal
import scipy.ndimage


def find_periodic_peaks(
    profile,
    refine=True,
    nfft=2**14,
    smooth_offset=4,
    num_offset_points=200,
    diagnostics=None,
):
    freqs, spectrum = scipy.signal.periodogram(
        profile, window="hann", nfft=nfft, scaling="spectrum"
    )
    pitch_idx = spectrum.argmax()
    pitch = 1 / freqs[pitch_idx]
    if diagnostics is not None:
        diagnostics["pitch"] = pitch
        spectrum_plot = hv.Curve(spectrum)
        spectrum_plot *= hv.VLine(pitch_idx).options(color="red")
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
        info = None
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
        info = pd.DataFrame({"hough_unshifted_value": profile[idxs], "shift": shifts})
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
    return refined_idxs, info
