#!/usr/bin/env python3

import os
import time
from tqdm import tqdm
import numpy
from scipy.signal import find_peaks
from skimage.filters import unsharp_mask
import tables
from multiprocessing import Pool
from params import read_params_file

"""
Write a PyTables table to each HDF5 FOV file, containing the trench co-ordinates for that FOV (and all time frames therein).
"""

# TODO add code to detect trench rows!
# TODO what needs to change if the images are now stored with F-ordering?


def run_trench_analysis(
    in_file,
    out_dir,
    filename,
    fov,
    frame,
    z_level,
    params,
    trench_width,
    trench_length,
    min_distance,
):
    """
    Run basic trench analysis on an HDF5 file: detect trenches for each time frame,
    or just using the first time frame. Then, measure trench properties in all time frames.

    Store detected trenches in a PyTables array.

    Each array entry (row) has 4 values (min_row, min_col, max_row, max_col), and there is 1 row per region (trench).
    Regions must be rectangular, and must have the same dimension (for stacking purposes).
    """
    h5file = tables.open_file(in_file, mode="r")
    node_path = "/{}/FOV_{}/Frame_{}/Z_{}/{}".format(
        filename, fov, frame, z_level, params["channel"]
    )

    if params["crop"]:
        img = h5file.get_node(node_path)[
            params["crop"]["top"] : params["crop"]["bottom"],
            params["crop"]["left"] : params["crop"]["right"],
        ]
    else:
        img = h5file.get_node(node_path).read()

    h5file.close()

    # Write coordinates of the detected trenches to HDF5
    path_string = "/{}/FOV_{}/Frame_{}/Z_{}".format(filename, fov, frame, z_level)
    out_file_path = os.path.join(out_dir, "{}/regions.h5".format(path_string))
    regions_file = tables.open_file(out_dir, mode="w")

    # Find all the trench rows in the image
    # Need to pass in the crop params, because these values have to be taken into consideration when returning co-ordinates
    # w.r.t. the original, uncropped image (which is the one which remains stored on disk).
    rows = detect_trench_rows(img, params["crop"])

    # For each row, find all the trenches
    # TODO does it make sense to store the rows all together, or separately?
    # Probably separately, because it might be the case that the trenches in different rows have different dimensions.
    for i, row in enumerate(rows):
        regions = detect_trenches(img, trench_width, trench_length, min_distance)
        # print("Detected {} trenches in File {} FOV {} Frame {} Z {} Row {}".format(regions.shape[0], filename, fov, frame, z_level, i)) # DEBUG
        regions_file.create_array(
            path_string, "row_{}".format(i), obj=regions, createparents=True
        )

    regions_file.close()


def detect_trench_rows(img, crop_params):
    """
    Find rows of trenches.
    Input a single image, and parameters [TBD].
    Return a numpy array, each row containing 4 co-ordinates:
    min_row, min_col, max_row, max_col, defining the rectangular region corresponding to that row of trenches.
    TODO sanity checking on the result, to make sure that it's within the boundaries of the image?
    FIXME what are the necessary parameters?
    """
    # List of lists: each sub-list has 4 values, and there is 1 sub-list per detected row
    coords = []

    # Flatten in the x-dimension
    peaks_data = img.mean(axis=1)

    # more steps TBD

    # Numpy conversion to 2D list from 2D list works as long as every sub-list has identical dimensions
    return numpy.array(coords)


def detect_trenches(img, trench_width, trench_length, min_distance):
    """
    Input an image, and some parameters.
    Output an numpy array of rectangular trench co-ordinates.
    """
    # TODO make sure that this is right
    img_width = img.shape[3] - img.shape[1]

    # Detect trenches (usually in phase contrast images).
    # Return (x_min, x_max, y_min, y_max) co-ordinates, for each trench,
    # OR min_row, min_col, max_row, max_col

    # Run the unsharp mask to enhance edges & to normalize from 0.0 - 1.0 # Values that worked: radius=2, amount=40
    img = unsharp_mask(
        img, radius=params["unsharp_mask_radius"], amount=params["unsharp_mask_amount"]
    )

    # Flatten in the y-dimension
    peaks_data = img.mean(axis=0)

    # Use scipy's peak detection algorithm
    peaks, _ = find_peaks(peaks_data)

    # Remove non-trench peaks
    trench_peaks = remove_non_trenches(peaks_data, peaks, params["cutoff"])

    # Merge peaks which are close to each other
    merged_peaks = merge_peaks(trench_peaks, min_distance)

    # Starting, ending x-coordinates for each trench.
    # These can be used to define rectangular regions of interest.
    # FIXME converting trench width to half width? best way?
    # If it's always the same, then pass this in rather than the whole width?
    trench_half_width = trench_width // 2
    ranges = define_trenches(merged_peaks, trench_half_width)

    # TODO 2 cleanup steps:

    # 1. Filter out overlapping trenches
    # UNIMPLEMENTED!

    # 2. Remove out-of-bounds trenches
    # TODO: replace if statments with while() loops, in case there are multiple trenches to be removed? (if they are overlapping...)
    # Remove first or last trench if they are out of bounds
    if len(ranges) != 0:
        if ranges[0][0] < 0:
            ranges.pop(0)

    # FIXME: input true image width, not this constant value
    if len(ranges) != 0:
        if ranges[-1][1] > img_width:
            ranges.pop()

    # Convert to numpy array
    return numpy.array(ranges)


def remove_non_trenches(peaks_data, peaks, cutoff):
    """
    Remove non-trenches.
    """
    # Filter out: only keep inter-trench peak indices
    threshold = (peaks_data[peaks] > cutoff) * peaks

    # Remove zero indices
    nz = numpy.nonzero(threshold)

    return peaks[nz]


def merge_peaks(peaks, min_peak_distance):
    """
    Merge the nearby peaks in peaks, a numpy array resulting from find_peaks().
    TODO: also add a sanity check: only merge them if both are above or below the threshold line
    """
    merged_peaks = []

    # Use double while loops to allow for merging clusters of nearby peaks
    i = 0
    while i < len(peaks) - 1:
        cluster = []
        cluster.append(peaks[i])

        j = i + 1
        # NOTE: comparison is always made to the i'th peak. This might behave weirdly
        # if there are many peaks, in series, which are spaced ~= min_peak_distance from one another
        while (j < len(peaks)) and (peaks[j] - peaks[i] < min_peak_distance):
            cluster.append(peaks[j])
            j += 1

        avg = numpy.mean(cluster).astype(numpy.int)
        merged_peaks.append(avg)
        i = j

    # Finally, the last peak
    # If i > len(peaks) - 1, then we merged the last one, and so can skip it,
    # but if i == len(peaks) - 1, then we have to keep it.
    if i == len(peaks) - 1:
        merged_peaks.append(peaks[i])

    # Convert from numpy list to numpy array
    return numpy.array(merged_peaks)


def define_trenches(trench_peaks, trench_half_width):
    """
    Detect Trenches
    Using just the trench midpoints & the pre-specified widths, define left -> right boundaries for each trench.
    Assume each trench spans from top to bottom.

    TODO: sanity bounds check near the edges of the FOV (when adding or subbing trench_half_width)
    Before making ranges, check whether would overlap?

    TODO: return 4 co-ordinates: either
    min_row, min_col, max_row, max_col,
    OR x_min, x_max, y_min, y_max
    """
    # ranges = []
    # for t in trench_peaks:
    # ranges.append( (t - trench_half_width, t + trench_half_width) )

    ranges = [(t - trench_half_width, t + trench_half_width) for t in trench_peaks]

    return ranges


def shared_region_linking(out_dir):
    """
    If regions are shared between time frames, then make symlinks within the HDF5 file
    so that each time frame still appears to have its own region data.
    """
    pass


def main_detection_function(out_dir, in_file, num_cpu, params_file, share_regions):
    """
    Main detection function.
    """
    h5file = tables.open_file(in_file, mode="r")

    # Iterate all the files & get their metadata
    # Loop the nodes immediately under Images to get the file names
    file_names = [i._v_name for i in h5file.list_nodes("/Images")]

    # Generate dict. of FOVs to process, for each ND2 file
    file_to_frames = {}

    # Progress bar
    total = 0

    # Metadata from each ND2 file
    metadata = {}
    if share_regions:
        for f in file_names:
            # Get metadata for each file
            n = h5file.get_node("/Metadata/{}".format(f))()
            m = get_metadata(n)

            # Total # of frames to be processed, for progress bar
            total += m["fields_of_view"] * m["z_levels"] * m["channels"]

            # Use the lowest numbered time frame to detect trenches in all time frames
            # (share the region across time frames).
            # NOTE is the sorting necessary, or do they always come out sorted from the conversion program?
            # Make it a list, so that the code in the main loop still works
            file_to_frames[f] = [sorted(m["frames"][0])]

            # store in a dictionary
            metadata[f] = m
    else:
        for f in file_names:
            # Get metadata for each file
            n = h5file.get_node("/Metadata/{}".format(f))()
            m = get_metadata(n)

            # Total # of frames to be processed, for progress bar
            total += m["fields_of_view"] * m["frames"] * m["z_levels"] * m["channels"]

            # Each time frame has its own detection (do not share the regions across time frames).
            # NOTE is the sorting necessary, or do they always come out sorted from the conversion program?
            file_to_frames[f] = sorted(m["frames"])

            # store in a dictionary
            metadata[f] = m

    h5file.close()

    # Read the parameters, and then convert the distances using pixel microns
    params = read_params_file(params_file)

    # Use the pixel microns for each file to convert distances to pixels
    file_trench_width = {}
    file_trench_length = {}
    file_min_distance = {}
    for f in file_names:
        pixel_microns = metadata[f]["pixel_microns"]
        # Dimensions (widths & lengths) of each trench
        file_trench_width[f] = int(params["trench_width"] / pixel_microns)
        file_trench_length[f] = int(params["trench_length"] / pixel_microns)
        # Minimum distance between trenches
        file_min_distance[f] = int(params["min_distance"] / pixel_microns)

    pbar = tqdm(total=total, desc="Frame #")

    def update_pbar(*a):
        pbar.update()

    # Run in parallel. Each File & FOV or Frame, can be processed independently.
    with Pool(processes=num_cpu) as p:
        for f in file_names:
            for fov in metadata[f]["fields_of_view"]:
                for frame in file_to_frames[f]:
                    # FIXME which z_level to use for trench detection? Probably makes sense to always share the same z_level,
                    # but that doesn't mean that the zeroth one should necessarily be the default.
                    # Could make this an additional parameter, and if None, then use the lowest one available.
                    z_level = 0  # (is there always a Z_0?)
                    args = [
                        in_file,
                        out_dir,
                        f,
                        fov,
                        frame,
                        z_level,
                        params,
                        file_trench_width[f],
                        file_trench_length[f],
                        file_min_distance[f],
                    ]
                    p.apply_async(run_trench_analysis, args, callback=update_pbar)
                    # run_trench_analysis(*args) # DEBUG

        pool.join()
        pool.close()

    # Done looping all files
    pbar.close()

    # TODO
    # FIXME What about linking into a super file, even if there aren't shared regions?
    # Linking, if share_regions
    if share_regions:
        shared_region_linking(out_dir)
