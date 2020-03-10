#!/usr/bin/env python3

import os
import time

import numpy

from scipy import ndimage  # For Gaussian blur in numpy / scipy

import mahotas
import skimage.morphology
import skimage.measure

import tables
import argparse
from multiprocessing import Pool


### Measure properties of trenches
### Update: include optional simple thresholding.
###
###

# Open the ND2 metadata table & extract the channel names (as a list)
def get_channel_names(in_dir, parent_file):
    h5file = tables.open_file(os.path.join(in_dir, parent_file), mode="r")
    metadata_table = h5file.get_node(
        "/tables/nd2_metadata"
    )()  # Last pair of parentheses necessary to dereference the node reference
    data = metadata_table.read()
    metadata_table.close()
    h5file.close()

    # Assume that there is only a single entry in this table
    # Returns bytes() (i.e. b'...'), which must be decoded() into a string
    channel_names = data[0]["info_channels"].decode().split("\t")
    return channel_names


# Define the column types for the PyTable: this stores un-segmented, whole trench information
def make_trench_type(channel_names):
    # Define the PyTables column types using a dictionary
    column_types = {
        "info_fov": tables.UInt16Col(),
        "info_frame": tables.UInt16Col(),
        "info_trench_number": tables.UInt16Col(),
    }

    # Add some column names dynamically
    # NOTE if using flat field correction, make these floats
    # Otherwise, integers
    for c in channel_names:
        column_types["total_intensity_{}".format(c)] = tables.UInt32Col()

    return type("Trench", (tables.IsDescription,), column_types)


# FIXME read ranges from the detected trenches table, extract from it the x_y dimensions & the stack size!!
def run_trench_analysis(
    fov_file_path,
    parent_file_path,
    identifier,
    channel_names,
    crop_top,
    crop_bottom,
    x_dimension,
    y_dimension,
    background,
):

    Trench = make_trench_type(channel_names)

    # Open the parent file to get the timestamps table
    parent_h5file = tables.open_file(parent_file_path, mode="r")
    table_timestamps = parent_h5file.get_node("/tables/fov_metadata")()

    # Open this FOV file
    fov_h5file = tables.open_file(fov_file_path, mode="r+")

    # Extract the FOV number from the file name
    base_name = os.path.basename(fov_file_path)
    (root, extension) = os.path.splitext(base_name)
    fov = int(root.split("_")[1])

    table_name = "measurements_whole_trench_{}".format(identifier)

    # Create the properties table
    # Create properties tables for each FOV, then merge them together at the end
    if fov_h5file.__contains__("/tables/{}".format(table_name)):
        node = fov_h5file.get_node("/tables/{}".format(table_name))
        node._f_remove()

    trench_properties_table = fov_h5file.create_table(
        "/tables", table_name, Trench, "Measurements_whole_trench"
    )
    trench_properties_row = trench_properties_table.row

    # Open the detected trenches table
    trenches_table_name = "/tables/trench_coordinates_{}".format(identifier)
    trenches_table = fov_h5file.get_node(trenches_table_name)

    ranges = [
        (row["bounding_box_min_col"], row["bounding_box_max_col"])
        for row in trenches_table.read()
    ]
    stack_size = len(ranges)

    # Loop all time points, re-use same trench locations
    for time_node in fov_h5file.root.images._f_iter_nodes():
        frame_number = int(time_node._v_name.split("_")[1])
        print("\t{}".format(time_node._v_name))

        # Only look at the zeroth returned match (there should only be 1 match)
        # NOTE could possibly speed up queries:
        # 1. indexing,
        # 2. or sort the iteration through nodes, and grab specific coordinates,
        # 3. or sort the iteration through nodes, and sort the columns
        timestamp_row = table_timestamps.read_where(
            """(info_fov == fov) & (info_frame == frame_number)"""
        )
        timestamp = timestamp_row["info_timestamp"][0]

        # First, load each of the channels & make trench stacks
        # Do this b/c need to refer across fluor. channels when making intensity measurements
        # FIXME actually, that's not necessary anymore? Does it matter?
        channel_stacks = {}
        for c in channel_names:
            channel_stacks[c] = load_stack_trenches(
                time_node._f_get_child(c),
                crop_top,
                crop_bottom,
                x_dimension,
                y_dimension,
                stack_size,
                ranges,
            )

        ### No segmentation
        # FIXME but who's the background for making a mask? Then we need to loop segmentation
        # channels all the same... or determine which is the appropriate channel empirically
        # and record which channel was used...
        # and then we also should store the masks!
        for n in range(channel_stacks[c].shape[2]):
            for c in channel_names:
                mask = channel_stacks[c][:, :, n] > background
                corrected_intensity = subtract_background_from_region(
                    channel_stacks[c][:, :, n], image_mask, background
                )
                trench_properties_row[
                    "total_intensity_{}".format(c)
                ] = corrected_intensity

                # trench_properties_row['total_intensity_{}'.format(c)] = channel_stacks[c][:, :, n].sum()

            trench_properties_row["info_fov"] = fov
            trench_properties_row["info_frame"] = frame_number
            trench_properties_row["info_trench_number"] = n

            # Append the trench information to the table
            trench_properties_row.append()

    ### Done measuring trenches in this FOV, across all time points
    table_timestamps.close()

    trench_properties_table.flush()
    trench_properties_table.close()

    trenches_table.close()

    fov_h5file.close()


# Load stack of image slices with dimensions pre-determined by detected trench positions.
# Also apply cropping to top, bottom of image, and apply constant background subtraction.
def load_stack_trenches(
    image_node,
    crop_top,
    crop_bottom,
    x_dimension,
    y_dimension,
    stack_size,
    ranges,
    background,
):
    stack = numpy.empty(shape=(y_dimension, x_dimension, stack_size), dtype=numpy.int16)
    # NOTE: Pass in the node reference, not the whole image, into this sub-routine.
    # Then, when using the slice notation below, will take advantage of the chunking so as to not ever
    # load or download the entire image.
    for i, r in enumerate(ranges):
        stack[:, :, i] = image_node[crop_top:crop_bottom, r[0] : r[1]]

    # FIXME can we use map() instead of a loop here?

    # Then, can apply operations to just the stack of trench cutouts
    # Background subtraction
    # Sets pixels > bg to their original values, < bg to bg, then remove bg.
    # stack = ( ((stack > background) * stack) + ((stack <= background) * background) ) - background
    # NOTE: this is not a good way to handle background. Rather, It's better to subtract it from
    # the final summed pixel intensities.

    return stack


# Input a list of pixel co-ordinates (as per skimage.measure.regionpprops),
# an associated image (2d numpy array) to read intensities from, and and a background signal value.
# Return the sum of the pixel intensities minus the product of the background value and the number of pixels.
def subtract_background_from_coords(coords, image, background):
    transposed = cell.coords.T
    summed_intensity = image[transposed[0], transposed[1]].sum()
    total_background = len(transposed) * background

    return summed_intensity - total_background


# Input an image, a binary mask, and a background value per pixel
# Return the sum of the pixel intensities minus the product of the background value and the number of pixels.
def subtract_background_from_region(image, image_mask, background):
    # Use mask to zero out unwanted pixels, and sum the intensities of the remaining ones
    summed_intensity = (image * image_mask).sum()

    # Calculate background signal
    # Sum of a mask is equivalent to the number of valid pixels, if the mask is as integer and not bool.
    total_background = image_mask.astype(numpy.uint8).sum() * background

    return summed_intensity - total_background


# TODO: throw errors if the arguments don't make sense?
def parse_args():
    parser = argparse.ArgumentParser(
        description="Measure fluorescence values within trenches, without cell segmentation."
    )
    parser.add_argument("-i", "--in-dir", type=str, help="Input HDF5 directory path")
    parser.add_argument(
        "-p", "--parent", type=str, help="Parent HDF5 file (with links & merged tables)"
    )
    parser.add_argument("-n", "--num-cpu", type=int, help="Number of CPUs to use")
    parser.add_argument("-I", "--identifier", type=str, help="Trench row identifier")
    parser.add_argument("-T", "--crop-top", type=int, help="Crop top")
    parser.add_argument("-B", "--crop-bottom", type=int, help="Crop bottom")
    parser.add_argument(
        "-t",
        "--threshold-file",
        type=int,
        help="File with pixel value minimum thresholds per channek.",
    )

    return parser.parse_args()


###

if __name__ == "__main__":
    ### Parse command-line arguments
    args = parse_args()

    in_dir = args.in_dir
    parent_file = args.parent
    num_cpu = args.num_cpu
    identifier = args.identifier
    crop_top = args.crop_top
    crop_bottom = args.crop_bottom
    threshold_file = args.threshold_file

    # FIXME what's the right way to read these in?
    x_dimension = 30
    y_dimension = crop_bottom - crop_top

    ###

    # Load the ND2 metadata to get the channel names
    channel_names = get_channel_names(in_dir, parent_file)

    # Load the thresholds from a text file
    # TODO

    ### Iterate all the files
    # NOTE: this could also be done outside of Python (e.g. in a Slurm / Bash script)
    # in an embarrasingly parallel fashion.

    print("Measuring trench properties...")
    start = time.time()

    fov_dir = os.path.join(in_dir, "FOV")
    files = os.listdir(fov_dir)
    ## DEBUG
    # files = ['FOV_0.h5']

    # Run in parallel
    # TODO: is it possible to pass in arguments non-explicitly, and therefore not need to use starmap?
    # Tradeoffs?
    if num_cpu > 1:
        args = [
            (
                os.path.join(fov_dir, filename),
                os.path.join(in_dir, parent_file),
                identifier,
                channel_names,
                crop_top,
                crop_bottom,
                x_dimension,
                y_dimension,
                threshold,
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
                    os.path.join(in_dir, parent_file),
                    identifier,
                    channel_names,
                    crop_top,
                    crop_bottom,
                    x_dimension,
                    y_dimension,
                    threshold,
                )

    ### Done looping all FOV
    end = time.time()
    print(end - start)
    print("Done measuring trench properties.")

    ### Merge the tables
    Trench = make_trench_type(channel_names)

    h5file = tables.open_file(os.path.join(in_dir, parent_file), mode="r+")
    table_name = "measurements_whole_trench_{}".format(identifier)

    # Create a table for storing the properties of the entire trench
    # If table exists, delete it and start over
    if h5file.__contains__("/tables/{}".format(table_name)):
        node = h5file.get_node("/tables/{}".format(table_name))
        node._f_remove()

    trench_properties_table_global = h5file.create_table(
        h5file.root.tables, table_name, Trench, "Measurements_whole_trench"
    )

    ### See: https://stackoverflow.com/questions/24891014/append-all-rows-from-one-table-to-another-using-pytables

    print("Merging tables...")
    start = time.time()

    for filename in os.listdir(fov_dir):
        if filename.endswith("h5"):
            table_file = tables.open_file(os.path.join(fov_dir, filename), "r")
            table = table_file.get_node("/tables/{}".format(table_name))

            for next_row in table.iterrows():
                trench_properties_table_global.append([next_row[:]])

            table.close()
            table_file.close()

    # Done copying trench properties
    trench_properties_table_global.flush()
    trench_properties_table_global.close()
    h5file.close()

    end = time.time()
    print(end - start)
    print("Done merging tables.")
