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


### Write a PyTables table to each HDF5 FOV file, containing the trench co-ordinates for that FOV
### (and all time frames therein)

# Run basic trench analysis on an HDF5 file:
# detect trenches using the first time frame, then measure trench properties in all time frames
# Store detected trenches in a PyTables  table.
# TODO also allow cropping left & right? (e.g. maybe it's blurry or something)
def run_trench_analysis(
    hdf5_file_path,
    min_peak_distance,
    cutoff,
    trench_half_width,
    crop_top,
    crop_bottom,
    identifier,
    channel_name,
    x_dimension,
    y_dimension,
):

    fov_h5file = tables.open_file(hdf5_file_path, mode="r+")

    ### Detect trenches in phase; return pairs of x_min, x_max co-ordinates, for each trench
    ranges = run_trench_detection(
        fov_h5file.root.images.Frame_0,
        channel_name,
        crop_top,
        crop_bottom,
        cutoff,
        trench_half_width,
    )
    # TODO: this will return the ndarray

    # FIXME instead of writing to a table, just use a ndarray with 4 uint16 values: x_min, x_max, y_min, y_max

    ### Create distinct properties tables for each FOV, then merge them together at the end
    # Coordinates of the detected trenches
    table_name = "trench_coordinates_{}".format(identifier)
    table_path = "/tables/{}".format(table_name)

    # TODO write a fixed-size array, no chunking.

    print("Detected {} trenches in FOV {}".format(ranges.shape[0], fov))


# node is the first frame in this fov, which will be used for trench detection
#
# Return pairs of x_min, x_max co-ordinates, for each trench
# Trenches must have identical dimensions
# FIXME pass in the true image width, not a constant value of 2048
def run_trench_detection(
    node, channel_name, crop_top, crop_bottom, cutoff, trench_half_width
):
    # TODO Write top left, bottom right corners to HDF5 table (for future reference) = 4 integers
    # TODO auto-detect cropping top & bottom

    # Get the very first phase image in this FOV position (zeroth time frame)
    img_phase_zeroth = node._f_get_child(channel_name)[crop_top:crop_bottom]

    # Run the unsharp mask to enhance edges & to normalize from 0.0 - 1.0
    # TODO tweak the radius & amount parameters? Pass in as variables?
    img_phase_zeroth = unsharp_mask(img_phase_zeroth, radius=2, amount=40)

    # Flatten in the y-dimension
    peaks_data = img_phase_zeroth.mean(axis=0)

    # Use scipy's peak detection algorithm
    peaks, _ = find_peaks(peaks_data)

    # Remove non-trench peaks
    trench_peaks = remove_non_trenches(peaks_data, peaks, cutoff)

    # Merge peaks which are close to each other
    merged_peaks = merge_peaks(trench_peaks, min_peak_distance)

    # Starting, ending x-coordinates for each trench.
    # These can be used to define rectangular regions of interest.
    ranges = define_trenches(merged_peaks, trench_half_width)

    ### 2 cleanup steps:

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
        if ranges[-1][1] > 2048:
            ranges.pop()

    return ranges


# Remove non-trenches
def remove_non_trenches(peaks_data, peaks, cutoff):
    # Filter out: only keep inter-trench peak indices
    threshold = (peaks_data[peaks] > cutoff) * peaks

    # Remove zero indices
    nz = numpy.nonzero(threshold)

    return peaks[nz]


# Merge the nearby peaks in peaks, a numpy array resulting from find_peaks()
# TODO: also add a sanity check: only merge them if both are above or below the threshold line
def merge_peaks(peaks, min_peak_distance):
    merged_peaks = []

    # Alternative version: double while loops
    # This allows for merging clusters of nearby peaks
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

    ## Done looping

    # Finally, the last peak
    # If i > len(peaks), then we merged the last one, and so can skip it
    # but if it's ==, then we have to keep it.
    if i == len(peaks) - 1:
        merged_peaks.append(peaks[i])

    ## Done with last peak

    # Convert to numpy array
    merged_peaks = numpy.array(merged_peaks)

    return merged_peaks


# Detect Trenches
def define_trenches(trench_peaks, trench_half_width):
    # Using just the trench midpoints & the pre-specified widths, define left-> right boundaries for each trench
    # Assume each trench spans from top to bottom

    # TODO: sanity bounds check near the edges of the FOV (when adding or subbing trench_half_width)
    # Before making ranges, check whether would overlap?

    # ranges = []
    # for t in trench_peaks:
    # ranges.append( (t - trench_half_width, t + trench_half_width) )

    ranges = [(t - trench_half_width, t + trench_half_width) for t in trench_peaks]

    return ranges


# TODO: instead of having these as separate tables, instead can *link* all of the output
# into a new HDF5 file, like we've been doing elsewhere.
# Merge the trench co-ordinates tables, store the result in the parent HDF5 file
def merge_tables(parent_file, in_dir, identifier):
    h5file = tables.open_file(os.path.join(in_dir, parent_file), mode="r+")

    # Create the table type
    Trench_Coords = make_trench_info_type()

    table_path = "/tables/trench_coordinates_{}".format(identifier)

    # Create a table for storing the properties of the entire trench
    # If table exists, delete it and start over
    if h5file.__contains__(table_path):
        node = h5file.get_node(table_path)
        node._f_remove()

    trench_coordinates_table_global = h5file.create_table(
        "/tables",
        "trench_coordinates_{}".format(identifier),
        Trench_Coords,
        "Trench_coordinates",
        createparents=True,
    )

    # For appending new entries
    trench_coordinates_row_global = trench_coordinates_table_global.row

    ### See: https://stackoverflow.com/questions/24891014/append-all-rows-from-one-table-to-another-using-pytables

    print("Merging tables...")
    start = time.time()

    directory = os.path.join(in_dir, "FOV")
    for filename in os.listdir(directory):
        if filename.endswith("h5"):
            table_file = tables.open_file(os.path.join(directory, filename), "r")
            local_table_path = "/tables/trench_coordinates_{}".format(identifier)
            table = table_file.get_node(table_path)

            for next_row in table.iterrows():
                trench_coordinates_table_global.append([next_row[:]])

            table.close()
            table_file.close()

    # Done copying trench properties
    trench_coordinates_table_global.close()
    h5file.close()

    ### Done merging!
    print("Done merging tables.")


### Main


def main_detection_function(out_dir, in_file, num_cpu, params_file, share_regions):
    params = read_params_file(params_file)

    ### Convert the parameters

    h5file_in = tables.open_file(in_file, mode="r")
    # Assume that all files have the same pixel_microns value
    # Pick the first available file, read its metadata, extract pixel_microns
    # FIXME
    node_list = []  # = h5file_in. ??
    pixel_microns = 1.0  # ? # TODO: read from file

    # FIXME: read pixel microns from the H5 file & convert the micron values to integer pixel numbers
    # Dimensions (widths & lengths) of each trench
    trench_width = int(params["trench_width"] / pixel_microns)
    trench_length = int(params["trench_length"] / pixel_microns)
    min_distance = int(params["min_distance"] / pixel_microns)

    ### Iterate all the files

    print("Detecting trenches...")
    start = time.time()

    # Run in parallel
    # Each File & FOV or Frame, can be processed independently.
    with Pool(processes=num_cpu) as p:
        if share_regions:
            # Files
            for filename in []:  # FIXME
                # FOVs
                for fov in h5file.iter_nodes("/"):
                    # Cannot assume that there is a Frame_0,
                    # but *can* sort the list of frames and pick out the lowest (earliest) one.
                    args = [params, frame_node, filename, fov, out_dir]
                    p.apply_async(run_trench_analysis, args)
                    # run_trench_analysis(*args) # DEBUG
        else:
            pass  # TODO: also iterate frames

    ### Done looping all FOV files
    end = time.time()
    print(end - start)
    print("Done detecting trenches.")

    # Merge all the detected regions
    merge_tables(parent_file, in_dir, identifier)
    # Done merging tables

    # FIXME there is still an HDF5 file not closed?
