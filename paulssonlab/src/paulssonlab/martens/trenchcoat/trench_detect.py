#!/usr/bin/env python3

import os
import time

import numpy

from scipy.signal import find_peaks
from skimage.filters import unsharp_mask

import tables
import argparse
from multiprocessing import Pool


### Write a PyTables table to each HDF5 FOV file, containing the trench co-ordinates for that FOV
### (and all time frames therein)
###


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

    # Create the table type
    # Have to do inside the sub-processes
    Trench_Coords = make_trench_info_type()

    ### Detect trenches in phase; return pairs of x_min, x_max co-ordinates, for each trench
    ranges = run_trench_detection(
        fov_h5file.root.images.Frame_0,
        channel_name,
        crop_top,
        crop_bottom,
        cutoff,
        trench_half_width,
    )

    ### Create distinct properties tables for each FOV, then merge them together at the end
    # Coordinates of the detected trenches
    table_name = "trench_coordinates_{}".format(identifier)
    table_path = "/tables/{}".format(table_name)

    # Remove the table if it already exists
    if fov_h5file.__contains__(table_path):
        node = fov_h5file.get_node(table_path)
        node._f_remove()

    # Create the table
    trench_coordinates_table = fov_h5file.create_table(
        "/tables", table_name, Trench_Coords, "Trench_coordinates", createparents=True
    )
    trench_coordinates_row = trench_coordinates_table.row

    # Extract the FOV number from the file name
    base_name = os.path.basename(hdf5_file_path)
    (root, extension) = os.path.splitext(base_name)
    fov = int(root.split("_")[1])

    # Write the trench co-ordinates into the table
    # TODO: determine if it's worth storing min_row, max_row... always the same?
    # TODO: use the cropped coordinates, or the original coordinates?
    for i, r in enumerate(ranges):
        # Write the trench coordinates to a table
        trench_coordinates_row["info_fov"] = fov
        trench_coordinates_row["info_trench_number"] = i
        # trench_coordinates_row['bounding_box_min_row'] = 0
        trench_coordinates_row["bounding_box_min_col"] = r[0]
        # trench_coordinates_row['bounding_box_max_row'] = crop_bottom - crop_top
        trench_coordinates_row["bounding_box_max_col"] = r[1]

        trench_coordinates_row.append()

    trench_coordinates_table.flush()
    trench_coordinates_table.close()

    print("Detected {} trenches in FOV {}".format(len(ranges), fov))


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


# 4 columns to store top left, bottom right corner coordinates of detected trenches,
# as well as FOV number & trench number.
def make_trench_info_type():
    # Define the PyTables column types using a dictionary
    column_types = {
        "info_fov": tables.UInt16Col(),
        "info_trench_number": tables.UInt16Col(),
        #         "bounding_box_min_row" : tables.UInt16Col(),
        "bounding_box_min_col": tables.UInt16Col(),
        #         "bounding_box_max_row" : tables.UInt16Col(),
        "bounding_box_max_col": tables.UInt16Col(),
    }

    return type("Trench_Coords", (tables.IsDescription,), column_types)


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


# TODO: throw errors if the arguments don't make sense?
def parse_args():
    parser = argparse.ArgumentParser(
        description="Detect trenches in phase contrast images."
    )
    parser.add_argument("-i", "--in-dir", type=str, help="Input HDF5 directory path")
    parser.add_argument(
        "-p", "--parent", type=str, help="Parent HDF5 file (with links & merged tables)"
    )
    parser.add_argument("-n", "--num-cpu", type=int, help="Number of CPUs to use")
    parser.add_argument(
        "-d",
        "--min-distance",
        type=int,
        help="Minimum peak distance (peaks which are at least this close will be merged)",
    )
    parser.add_argument(
        "-c", "--cutoff", type=float, help="Normalized trench cutoff (0 - 1)"
    )
    parser.add_argument(
        "-W", "--half-width", type=int, help="Half the width of each trench"
    )
    parser.add_argument("-T", "--crop-top", type=int, help="Crop top")
    parser.add_argument("-B", "--crop-bottom", type=int, help="Crop bottom")
    parser.add_argument("-I", "--identifier", type=str, help="Trench row identifier")
    parser.add_argument(
        "-C", "--channel", type=str, help="Channel to use for trench detection"
    )

    return parser.parse_args()


### Main

if __name__ == "__main__":
    ### Parse command-line arguments
    args = parse_args()

    # 1. HDF directory, which contains an FOV directory
    # Store a separate HDF5 file for each FOV in this directory
    in_dir = args.in_dir

    # 2. The parent HDF5 file
    parent_file = args.parent

    # 3. Number of processes to run in parallel (num-cpu)
    num_cpu = args.num_cpu

    # 4. Peaks which are at least this close will be merged.
    min_peak_distance = args.min_distance  # 25

    # 5. A peak less than this value is in between trenches,
    # greater than this value is a trench.
    cutoff = args.cutoff  # 0.35

    # 6. Manually specify trench widths
    trench_half_width = args.half_width  # 15

    # 7. Remove extraneous regions from the image
    # TODO: automatically determine cropped area
    crop_top = args.crop_top  # 1000
    crop_bottom = args.crop_bottom  # 1325

    # 8. Specify a unique identifier for this set of trenches.
    # Identifiers are useful if multiple rows are to be designated,
    # or if different detection algorithms are to be compared.
    # Example identifier could be a number, or a detection parameter
    identifier = args.identifier

    # 9. The name of the channel to use for trench detection
    channel_name = args.channel

    ###

    # Dimensions (widths & lengths) of each trench
    x_dimension = trench_half_width * 2
    y_dimension = crop_bottom - crop_top

    ## Instantiate the table row data types
    # segmentation_enum = tables.Enum(segmentation_channels)

    ### Iterate all the files

    print("Detecting trenches...")
    start = time.time()

    fov_dir = os.path.join(in_dir, "FOV")

    files = os.listdir(fov_dir)
    ## DEBUG
    # files = ["FOV_0.h5"]

    # Run in parallel
    if num_cpu > 1:
        args = [
            (
                os.path.join(fov_dir, filename),
                min_peak_distance,
                cutoff,
                trench_half_width,
                crop_top,
                crop_bottom,
                identifier,
                channel_name,
                x_dimension,
                y_dimension,
            )
            for filename in files
        ]

        # TODO: weed out any files which do not end in h5

        with Pool(processes=num_cpu) as p:
            p.starmap(run_trench_analysis, args)

    # Do not run in parallel
    else:
        for filename in files:
            if filename.endswith("h5"):
                run_trench_analysis(
                    os.path.join(fov_dir, filename),
                    min_peak_distance,
                    cutoff,
                    trench_half_width,
                    crop_top,
                    crop_bottom,
                    identifier,
                    channel_name,
                    x_dimension,
                    y_dimension,
                )

    ### Done looping all FOV files
    end = time.time()
    print(end - start)
    print("Done detecting trenches.")

    # TODO: should we then merge the detected trench tables?

    # FIXME input file names
    merge_tables(parent_file, in_dir, identifier)
    # Done merging tables

    # FIXME there is still an HDF5 file not closed?
