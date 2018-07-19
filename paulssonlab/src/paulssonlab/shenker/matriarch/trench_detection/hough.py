import numpy as np
import pandas as pd
import scipy
import skimage
import skimage.morphology
import skimage.segmentation
from itertools import zip_longest
import holoviews as hv
import holoviews.operation.datashader as datashader
# TODO: fix imports
from holoborodko_diff import holo_diff
from .core import (edge_point, coords_along,
                   stack_jagged, stack_jagged_points)
from .periodogram import label_for_trenches
from util import getitem_if_not_none
from ui import RevImage
from image import (hough_line_intensity, remove_large_objects,
                   normalize_componentwise, gaussian_box_approximation)
from geometry import get_image_limits
from workflow import points_dataframe
import common

def find_periodic_lines(img, theta=None, smooth=4, hough_func=hough_line_intensity,
                        num_offset_points=200,
                        interpolate=5,
                        refine=None, diagnostics=None):
    if theta is None:
        theta = np.linspace(np.deg2rad(-45), np.deg2rad(45), 90)
    h, theta, rho = hough_func(img, theta=theta)
    diff_h = np.diff(h.astype(np.int_), axis=1) # TODO: is diff necessary??
    diff_h_std = diff_h.std(axis=0)# / diff_h.max(axis=0)
    if smooth:
        diff_h_std_smoothed = scipy.ndimage.filters.gaussian_filter1d(diff_h_std, smooth)
    else:
        diff_h_std_smoothed = diff_h_std
    theta_idx = diff_h_std_smoothed.argmax()
    angle = theta[theta_idx]
    if diagnostics is not None:
        diagnostics['input'] = RevImage(img)
        diagnostics['angle_range'] = (np.rad2deg(theta[0]), np.rad2deg(theta[-1]))
        bounds = (np.rad2deg(theta[0]), rho[0], np.rad2deg(theta[-1]), rho[-1])
        diagnostics['log_hough'] = hv.Image(np.log(1+h), bounds=bounds)
        # TODO: fix the left-edge vs. right-edge issue for theta bins
        theta_degrees = np.rad2deg(theta[:-1])
        diff_h_plot = hv.Curve((theta_degrees, diff_h_std)) 
        if smooth:
            diff_h_plot *= hv.Curve((theta_degrees, diff_h_std_smoothed)).options(color='cyan')
        diff_h_plot *= hv.VLine(np.rad2deg(angle)).options(color='red')
        diagnostics['diff_h_std'] = diff_h_plot
        diagnostics['angle'] = np.rad2deg(angle)
    profile = h[:,theta_idx]
    # TODO: replace with trim by rhos
    x_lim, y_lim = get_image_limits(img.shape)
    x_min, x_max = x_lim
    y_min, y_max = y_lim
    ############
    max_dist = int(np.ceil(np.sqrt(x_max**2 + y_max**2)))
    #delta = 1024
    rho_min = #-y_max*np.sin(angle) + delta
    rho_max = #x_max*np.cos(angle) + delta
    print('rho range',rho_min,rho_max)
    ############
    trimmed_profile = np.trim_zeros(profile)
    #if diagnostics is not None:
    #    diagnostics['trimmed_profile'] = hv.Curve(trimmed_profile)
    if interpolate:
        interpolated_profile = trimmed_profile
    else:
        interpolated_profile = trimmed_profile
    freqs, spectrum = scipy.signal.periodogram(interpolated_profile, scaling='spectrum')
    #spectrum[:2] = 0 # TODO: ignore two lowest frequencies
    pitch_idx = spectrum.argmax()
    pitch = 1/freqs[pitch_idx]
    if diagnostics is not None:
        diagnostics['pitch'] = pitch
        diagnostics['spectrum'] = hv.Curve(spectrum) * hv.VLine(pitch_idx).options(color='red')
    offsets = np.linspace(0, pitch, num_offset_points, endpoint=False)
    # always leave a period of length `pitch` so we can add offsets
    # even when we could've fit another period
    offset_idxs = (np.arange(0, len(profile), pitch)[:-1] + offsets[:,np.newaxis]).astype(np.int_)
    #anchor_rhos = _anchor_rhos(angle, pitch, 0, x_lim, y_lim)
    #print('$',anchor_rhos)
    #anchor_idxs = np.where(anchor_rhos == rho[:,np.newaxis])[0]
    #print('>',anchor_idxs)
    #max_dist = int(np.ceil(np.sqrt(x_max**2 + y_max**2)))
    #offset_idxs = (max_dist + anchor_rhos + offsets[:,np.newaxis]).astype(np.int_)
    #print('>>',offset_idxs)
    offset_objective = profile[offset_idxs].sum(axis=1)
    if smooth:
        offset_objective_smoothed = scipy.ndimage.filters.gaussian_filter1d(offset_objective, smooth)
    else:
        offset_objective_smoothed = offset_objective
    offset_idx = offset_objective_smoothed.argmax()
    offset = offsets[offset_idx]
    if diagnostics is not None:
        diagnostics['offset'] = offset
        offset_plot = hv.Curve((offsets, offset_objective))
        if smooth:
            offset_plot *= hv.Curve((offsets, offset_objective_smoothed)).options(color='cyan')
        offset_plot *= hv.VLine(offset).options(color='red')
        diagnostics['offsets'] = offset_plot
    
        periodic_points = hv.Scatter((rho[offset_idxs[offset_idx]], profile[offset_idxs[offset_idx]])).options(size=5, color='red')
        profile_plot = hv.Curve((rho, profile)) * periodic_points
        profile_plot *= hv.VLine(rho_min).options(color='red')
        profile_plot *= hv.VLine(rho_max).options(color='red')
        if refine:
            refined_points = hv.Scatter((rho[offset_idxs[offset_idx]], profile[offset_idxs[offset_idx]])).options(size=5, color='cyan')
            profile_plot *= refined_points
        diagnostics['profile'] = profile_plot
    return angle, pitch, offset

def find_trench_lines(img, window=np.deg2rad(10), refine=0.5, diagnostics=None):
    angle1, pitch1, offset1 = find_periodic_lines(img,
                                                  refine=None,
                                                  diagnostics=getitem_if_not_none(diagnostics, 'hough_1'))
    angle2, pitch2, offset2 = find_periodic_lines(img,
                                                  theta=np.linspace(angle1-window, angle1+window, 200),
                                                  refine=refine,
                                                  diagnostics=getitem_if_not_none(diagnostics, 'hough_2'))
    return angle2, pitch2, offset2

def _anchor_rhos(angle, pitch, offset, x_lim, y_lim):
    x_min, x_max = x_lim
    y_min, y_max = y_lim
    max_dist = int(np.ceil(np.sqrt(x_max**2 + y_max**2)))
    # TODO: not totally clear that the signs are right on max_dist, but seems to work
    if angle < 0:
        effective_offset = max_dist + x_max*np.cos(angle) - offset
    else:
        effective_offset = -max_dist + offset
    abs_angle = np.abs(angle)
    delta = (y_max - x_max*np.tan(abs_angle))*np.sin(abs_angle)
    rhos = np.arange(effective_offset % pitch,
                     x_max/np.cos(angle) + delta, pitch)
    return rhos

def trench_anchors(angle, rhos, x_lim, y_lim):
    x_min, x_max = x_lim
    y_min, y_max = y_lim
    anchors = rhos[:,np.newaxis] * np.array((np.cos(angle), np.sin(angle)))[np.newaxis,:]
    if angle < 0:
        upper_right = np.array((x_max,0))
        anchors = upper_right - anchors
    return anchors

def _anchor_rhos2(angle, pitch, offset, x_lim, y_lim):
    x_min, x_max = x_lim
    y_min, y_max = y_lim
    max_dist = int(np.ceil(np.sqrt(x_max**2 + y_max**2)))
    # TODO: not totally clear that the signs are right on max_dist, but seems to work
    if angle < 0:
        effective_offset = max_dist + x_max*np.cos(angle) - offset
    else:
        effective_offset = -max_dist + offset
    abs_angle = np.abs(angle)
    delta = (y_max - x_max*np.tan(abs_angle))*np.sin(abs_angle)
    rhos = np.arange(effective_offset % pitch,
                     x_max/np.cos(angle) + delta, pitch)
    return rhos

def trench_anchors2(angle, pitch, offset, x_lim, y_lim):
    x_min, x_max = x_lim
    y_min, y_max = y_lim
    rhos = _anchor_rhos2(angle, pitch, offset, x_lim, y_lim)
    anchors = rhos[:,np.newaxis] * np.array((np.cos(angle), np.sin(angle)))[np.newaxis,:]
    if angle < 0:
        upper_right = np.array((x_max,0))
        anchors = upper_right - anchors
    return anchors

def find_trench_ends(img, angle, pitch, offset, margin=15,
                     threshold=0.8, smooth=100,
                     diagnostics=None):
    x_lim, y_lim = get_image_limits(img.shape)
    anchors = trench_anchors(angle, pitch, offset, x_lim, y_lim)
    profiles = []
    line_points = []
    offsets = []
    for anchor in anchors:
        top_anchor = edge_point(anchor, 3/2*np.pi-angle, x_lim, y_lim)
        bottom_anchor = edge_point(anchor, np.pi/2-angle, x_lim, y_lim)
        line_length = np.linalg.norm(top_anchor-bottom_anchor)
        top_length = np.linalg.norm(top_anchor-anchor)
        bottom_length = np.linalg.norm(bottom_anchor-anchor)
        xs, ys = coords_along(top_anchor, bottom_anchor)
        profile = img[ys,xs]
        points = np.vstack((xs, ys)).T
        # TODO: precision??
        if line_length >= max(top_length, bottom_length):
            # line contains anchor
            offset = -int(np.ceil(top_length))
        else:
            # line is strictly on one side of anchor
            if top_length <= bottom_length:
                # line lies below anchor
                offset = int(np.ceil(top_length))
            else:
                # line lies above anchor
                offset = int(np.ceil(bottom_length))
        profiles.append(profile)
        line_points.append(points)
        offsets.append(offset)
    min_offset = min(offsets)
    anchor_idx = -min_offset
    max_stacked_length = max([len(profile)-offset for offset, profile in zip(offsets, profiles)])
    padded_profiles = []
    padded_line_points = []
    for i, (profile, points, offset) in enumerate(zip(profiles, line_points, offsets)):
        left_padding = offset - min_offset
        right_padding = max_stacked_length - left_padding - len(profile)
        padded_profile = np.pad(profile, (left_padding, right_padding), 'constant', constant_values=np.nan)
        padded_points = np.pad(points, [(left_padding, right_padding), (0,0)], 'edge')
        padded_profiles.append(padded_profile)
        padded_line_points.append(padded_points)
    if diagnostics is not None:
        # TODO: make hv.Path??
        diagnostics['profiles'] = hv.Overlay.from_values([hv.Curve(tp) for tp in padded_profiles])
    #stacked_profile = padded_profiles[25]
    stacked_profile = np.nanpercentile(np.array(padded_profiles), threshold*100, axis=0)
    stacked_points = np.array(padded_line_points).swapaxes(0,1)
    if diagnostics is not None:
        diagnostics['threshold'] = threshold
        lines_plot = hv.Path([[points[0], points[-1]] for points in line_points]).options(color='blue')
        top_line_plot = hv.Points([points[0] for points in line_points]).options(color='green')
        bottom_line_plot = hv.Points([points[-1] for points in line_points]).options(color='red')
        anchor_points_plot = hv.Points(anchors).options(size=3, color='cyan')
        diagnostics['image_with_lines'] = RevImage(img) * \
            lines_plot * \
            top_line_plot * \
            bottom_line_plot * \
            anchor_points_plot
    stacked_profile_diff = holo_diff(1, stacked_profile)
    # using np.nanargmax/min because we might have an all-nan axis
    top_end = max(np.nanargmax(stacked_profile_diff) - margin, 0)
    bottom_end = min(np.nanargmin(stacked_profile_diff) + margin, len(stacked_profile)-1)
    if diagnostics is not None:
        diagnostics['margin'] = margin
        diagnostics['stacked_profile'] = hv.Curve(stacked_profile) * \
            hv.Curve(stacked_profile_diff).options(color='cyan') * \
            hv.VLine(anchor_idx).options(color='gray') * \
            hv.VLine(top_end).options(color='green') * \
            hv.VLine(bottom_end).options(color='red')
    top_endpoints = stacked_points[top_end]
    bottom_endpoints = stacked_points[bottom_end]
    # discard trenches where top endpoint is the same as the bottom endpoint
    mask = ~np.apply_along_axis(np.all, 1, np.equal(top_endpoints, bottom_endpoints))
    top_endpoints = top_endpoints[mask]
    bottom_endpoints = bottom_endpoints[mask]
    if diagnostics is not None:
        trench_plot = hv.Path([[top_endpoint, bottom_endpoint]
                                  for top_endpoint, bottom_endpoint
                                  in zip(top_endpoints, bottom_endpoints)]).options(color='white')
        top_points_plot = hv.Points(top_endpoints).options(size=3, color='green')
        bottom_points_plot = hv.Points(bottom_endpoints).options(size=3, color='red')
        diagnostics['image_with_trenches'] = RevImage(img) * \
            trench_plot * \
            top_points_plot * \
            bottom_points_plot
    df_columns = {'top': points_dataframe(top_endpoints),
                  'bottom': points_dataframe(bottom_endpoints)}
    df = pd.concat(df_columns, axis=1)
    return df

def find_trenches(img, setwise=True, diagnostics=None):
    img_normalized, img_labels, label_index = label_for_trenches(img, diagnostics=getitem_if_not_none(diagnostics, 'labeling'))
    trench_sets = {}
    if not setwise:
        img_masked = np.where(skimage.morphology.binary_dilation(img_labels != 0), img_normalized, np.percentile(img_normalized, 5))
        angle, pitch, offset = find_trench_lines(img_masked, diagnostics=getitem_if_not_none(diagnostics, 'find_trench_lines'))
    for label in label_index:
        label_diagnostics = getitem_if_not_none(diagnostics, 'label_{}'.format(label))
        img_masked = np.where(skimage.morphology.binary_dilation(img_labels == label), img_normalized, np.percentile(img_normalized, 5))
        if setwise:
            angle, pitch, offset = find_trench_lines(img_masked, diagnostics=getitem_if_not_none(label_diagnostics, 'find_trench_lines'))
        trench_sets[label] = find_trench_ends(img_masked, angle, pitch, offset, diagnostics=getitem_if_not_none(label_diagnostics, 'find_trench_ends'))
    trenches_df = pd.concat(trench_sets)
    trenches_df.index.names = ['trench_set', 'trench']
    return trenches_df