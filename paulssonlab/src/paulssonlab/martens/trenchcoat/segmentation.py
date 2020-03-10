#!/usr/bin/env python3

import os
import time

import numpy

from scipy import ndimage

import mahotas
import skimage.morphology
import skimage.measure

import tables
import argparse
from multiprocessing import Pool


### Measure fluorescence intensities in trenches after performing cell segmentation.
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


# Input path to parameters file, return dictionary of parameters
def read_params_file(file_path):
    # Input channel name, output list of segmentation parameters
    params = {}

    with open(file_path, "r") as filehandle:
        first_line = filehandle.readline().strip()
        param_names = first_line.split("\t")

        for line in filehandle:
            line = line.strip()
            elements = line.split("\t")

            # Get the channel name from the first element, then remove from the list
            channel_name = elements[0]
            del elements[0]

            # Input parameter name, output parameter value
            these_params = {}
            for i, e in enumerate(elements):
                # i + 1 b/c channel name was not removed from this list
                name = param_names[i + 1]
                # Some params are numeric...
                if name == "threshold":
                    if e != "-":
                        these_params[name] = int(e)
                elif name == "background":
                    if e != "-":
                        these_params[name] = int(e)
                elif name == "otsu_multiplier":
                    if e != "-":
                        these_params[name] = float(e)
                elif name == "niblack_k":
                    if e != "-":
                        these_params[name] = float(e)
                elif name == "niblack_w":
                    if e != "-":
                        these_params[name] = int(e)
                else:
                    these_params[name] = e

            params[channel_name] = these_params

    return params


# Input params dict, return list of fluorescence channels
def parse_fluorescence_channel_names(params):
    fluor_channels = []
    for key, value in params.items():
        if params[key]["BF_or_FLUOR"] == "FLUOR":
            fluor_channels.append(key)

    return fluor_channels


# Input params dict, return list of segmentation channels
def parse_segmentation_channel_names(params):
    segmentation_channels = []
    for key, value in params.items():
        if params[key]["is_segmentation"] == "True":
            segmentation_channels.append(key)

    return segmentation_channels


# Define the column types for the PyTable: this stores segmented cell information
def make_cell_type(channels, seg_channels):
    # Define the PyTables column types using a dictionary
    column_types = {
        "info_fov": tables.UInt16Col(),
        "info_frame": tables.UInt16Col(),
        "info_seg_channel": tables.EnumCol(seg_channels, seg_channels[0], "int8"),
        "info_label": tables.UInt16Col(),
        "info_trench_number": tables.UInt16Col(),
        "geometry_Area": tables.UInt16Col(),
        "geometry_Orientation": tables.Float32Col(),
        "geometry_Perimeter": tables.Float32Col(),
        # See note below about math domain errors.
        "axis_length_major": tables.Float32Col(),
        "axis_length_minor": tables.Float32Col(),
        "centroid_row": tables.UInt16Col(),
        "centroid_col": tables.UInt16Col(),
        "bounding_box_min_row": tables.UInt16Col(),
        "bounding_box_min_col": tables.UInt16Col(),
        "bounding_box_max_row": tables.UInt16Col(),
        "bounding_box_max_col": tables.UInt16Col(),
    }

    # Add some column names dynamically
    # NOTE if using flat field correction, make these floats
    # Otherwise, integers
    for c in channels:
        # Total intensities
        column_types["total_intensity_{}".format(c)] = tables.UInt32Col()

        # Intensity ratios
        column_types["ratio_{}_vs_seg".format(c)] = tables.Float32Col()

    return type("Cell", (tables.IsDescription,), column_types)


# Very basic thresholding!
def run_thresholding(stack, threshold):
    return (stack > threshold).astype(numpy.uint8)


# Thresholding which also uses Niblack (mean, std dev) & watershedding
def run_segmentation(stack, otsu_multiplier, niblack_k, niblack_w, min_size):
    scaling_factor = 10
    garbage_otsu_value = 110

    # Store a stack of labeled segmentation masks
    result = numpy.empty(stack.shape, dtype=numpy.int8)

    # For each image in the stack
    for z in range(stack.shape[2]):
        # Grab the img from the stack
        stack_elem = stack[:, :, z]

        # Threshold is calculated using Otsu method, which samples the distribution of pixel intensities, and then
        # scaled by the otsu multiplier, which might differ depending on empirical experience with a given channel or dataset.
        threshold = skimage.filters.threshold_otsu(stack_elem) * otsu_multiplier

        # If it's a low value, then it's probably just noise...
        if threshold < garbage_otsu_value:
            result[:, :, z] = numpy.zeros(shape=(stack.shape[0], stack.shape[1]))
        else:
            # Set pixels < threshold to zero, but preserve values > threshold.
            # Helps with the mean + stdev*k calculations in Niblack.
            # It is "bimodal" because the values are either zero, or they are the original value.
            bimodal = (stack_elem < threshold) * 0 + (
                stack_elem >= threshold
            ) * stack_elem.astype(numpy.uint32) * scaling_factor

            # The multiply by scaling_factor trick helps to emphasize the value of the pixels which are kept.
            # This enhances the std. dev. at the edges of cells.

            # Apply Niblack method
            niblack = skimage.filters.threshold_niblack(
                bimodal, window_size=niblack_w, k=niblack_k
            )

            # Only keep pixels which are < than the mean*stdev*k (in rectangular region, size w).
            mask = stack_elem < niblack

            # Closing helps to fill in internal holes inside the cell
            mask = skimage.morphology.binary_closing(mask)

            # Remove small objects
            skimage.morphology.remove_small_objects(
                mask, min_size=min_size, connectivity=2, in_place=True
            )

            # Label the regions
            markers = skimage.measure.label(mask)

            # Make a basin (invert max, min pixels)
            # NOTE it would be nice if the watershed algo. allowed for a max priority queue.
            image = stack_elem.max() - stack_elem

            result[:, :, z] = skimage.morphology.watershed(
                image=image, markers=markers, mask=mask
            )

    return result


# Load stack of image slices with dimensions pre-determined by detected trench positions.
# Also apply cropping to top, bottom of image
def load_stack_trenches(
    image_node, crop_top, crop_bottom, x_dimension, y_dimension, stack_size, ranges
):
    stack = numpy.empty(
        shape=(y_dimension, x_dimension, stack_size), dtype=numpy.uint16
    )
    # NOTE: Pass in the node reference, not the whole image, into this sub-routine.
    # Then, when using the slice notation below, will take advantage of the chunking so as to not ever
    # load or download the entire image.
    for i, r in enumerate(ranges):
        stack[:, :, i] = image_node[crop_top:crop_bottom, r[0] : r[1]]

    # FIXME can we use map() instead of a loop here? Best way to loop over z in numpy array?

    return stack


# Open the table of trench coordinates for this FOV & return the coordinates as a list of pairs (min_col, max_col).
# TODO: also include the min_row & max_row information, always?
def get_trench_coordinates(path, identifier, h5file):
    trenches_table_name = "{}_{}".format(path, identifier)
    trenches_table = h5file.get_node(trenches_table_name)
    ranges = [
        (row["bounding_box_min_col"], row["bounding_box_max_col"])
        for row in trenches_table.read()
    ]
    trenches_table.close()

    return ranges


# Delete the properties table if it already exists
# Initialize it
# Return the row accessor for appending new entries
# FIXME name vs title? What to pass in, and how?
def prepare_properties_table(table_name, Cell, h5file):
    path = "/tables/{}".format(table_name)

    if h5file.__contains__(path):
        node = h5file.get_node(path)
        node._f_remove()

    properties_table = h5file.create_table(
        "/tables", table_name, Cell, "Measurements_thresholded_trench"
    )
    properties_row = properties_table.row

    return (properties_table, properties_row)


# TODO: is there a way to merge some of this code with the whole-trench code? Reduce duplication...
# FIXME read ranges from the detected trenches table, extract from it the x_y dimensions & the stack size!!
# NOTE not true segmentation; use the terms segmentation & thresholding somewhat interchangeably...
def run_segmentation_analysis(
    fov_file_path,
    identifier,
    channels,
    seg_channels,
    crop_top,
    crop_bottom,
    x_dimension,
    y_dimension,
    params,
):

    # ### HDF5 Preparation

    # For writing the segmentation masks
    filters = tables.Filters(complevel=1, complib="zlib")
    chunk_dimensions = (
        32,
        32,
    )  # Smaller chunks probably work ok, since the trenches are small?
    # NOTE if segmenting trenches, is it even worth chunking? Does it matter?

    # TODO Is there a way to pass in the function, if we can't pass in the type?
    Cell = make_cell_type(channels, seg_channels)

    # Define the enumerable type, for storing in the table
    seg_enum = tables.Enum(seg_channels)

    # Open this FOV file
    fov_h5file = tables.open_file(fov_file_path, mode="r+")

    # Extract the FOV number from the file name
    base_name = os.path.basename(fov_file_path)
    (root, extension) = os.path.splitext(base_name)
    fov = int(root.split("_")[1])

    print("FOV_{}".format(fov))

    # Create the properties table
    # Create properties tables for each FOV, then merge them together at the end
    table_name = "measurements_thresholded_trench_{}".format(identifier)
    (properties_table, properties_row) = prepare_properties_table(
        table_name, Cell, fov_h5file
    )

    # And if the masks already exist, then delete them too
    if fov_h5file.__contains__("/masks_{}".format(identifier)):
        node = fov_h5file.get_node("/masks_{}".format(identifier))
        node._f_remove(recursive=True)

    # Open the detected trenches table
    ranges = get_trench_coordinates(
        "/tables/trench_coordinates", identifier, fov_h5file
    )
    stack_size = len(ranges)

    # ### Done preparing the HDF5 file & tables

    min_size = 40

    # Loop all time points, re-use same trench locations
    for time_node in fov_h5file.root.images._f_iter_nodes():
        frame_number = int(time_node._v_name.split("_")[1])

        ### Load a stack of trenches into a dict:
        # Input channel name, return stack of trenches
        channel_to_stack = {}
        for c in channels:
            channel_to_stack[c] = load_stack_trenches(
                time_node._f_get_child(c),
                crop_top,
                crop_bottom,
                x_dimension,
                y_dimension,
                stack_size,
                ranges,
            )

        ### Run segmentation & write mask & properties to disk
        for sc in seg_channels:

            ## Threshold  the seg. channel stack, create a mask stack
            # mask_stack = run_thresholding(channel_to_stack[sc], params[sc]['threshold'])

            mask_stack = run_segmentation(
                channel_to_stack[sc],
                params[sc]["otsu_multiplier"],
                params[sc]["niblack_k"],
                params[sc]["niblack_w"],
                min_size,
            )

            ## Iterate each trench in the stack, compute the properties, & write to a table
            for z in range(stack_size):
                ## Store the thresholding result (the "mask") in the HDF5 file
                # NOTE even empty thresholding masks will be stored. Is this necessary?
                fov_h5file.create_carray(
                    "/masks_{}/Frame_{}/{}".format(identifier, frame_number, sc),
                    "trench_{}".format(z),
                    obj=mask_stack[:, :, z],
                    title="{}/{} thresholding mask".format(sc, z),
                    chunkshape=chunk_dimensions,
                    filters=filters,
                    createparents=True,
                )

                # Compute the properties
                properties = skimage.measure.regionprops(mask_stack[:, :, z])

                # Get the enumerable type for table columns
                sc_to_enum = seg_enum[sc]

                # TODO: if there is an underflow error, how to handle it? (we're working with uint16)
                #       1. The best thing to do would be to find a way to write NaN. However, we're storing ints, not floats.
                #       2. Or add a column which flags if there was underflow.
                #       3. Or skip writing the trench information entirely.
                #          The absence of the datapoint itself would be the signal that there was a problem.
                #       4. Or store 32-bit signed ints, or floats, and allow for negative numbers.

                for region in properties:

                    write_properties_to_table(
                        fov,
                        frame_number,
                        region,
                        sc,
                        channel_to_stack,
                        params,
                        properties_row,
                        sc_to_enum,
                        z,
                        channels,
                    )

    ### Done measuring trenches in this FOV, across all time points
    properties_table.flush()
    properties_table.close()

    fov_h5file.close()

    ### Done!


# Once a given trench has a segmentation mask, then write the properties of the masked region to the table.
def write_properties_to_table(
    fov,
    frame_number,
    region,
    sc,
    channel_to_stack,
    params,
    properties_row,
    enum,
    z,
    channels,
):

    # NOTE skimage is changing their coordinates -- do we still want to transpose??? I think so...
    coords = region.coords.T

    # First, the seg. channel. Do it first b/c we want it for the ratio.
    # NOTE is it faster to use the coords[] array, or to "multiply" by the mask?
    corrected_intensity_seg = subtract_background_from_coords(
        coords, channel_to_stack[sc][:, :, z], params[sc]["background"]
    )
    properties_row["total_intensity_{}".format(sc)] = corrected_intensity_seg

    # Then, the other channels
    for c in channels:
        if c != sc:
            corrected_intensity = subtract_background_from_coords(
                coords, channel_to_stack[c][:, :, z], params[c]["background"]
            )

            properties_row["total_intensity_{}".format(c)] = corrected_intensity

            # And calculate a ratio:
            # TODO do we ever get NaN or divide by zero here? (If there are no pixels which meet the threshold...)
            if corrected_intensity_seg == 0:
                properties_row["ratio_{}_vs_seg".format(c)] = numpy.nan
            else:
                properties_row[
                    "ratio_{}_vs_seg".format(c)
                ] = corrected_intensity.astype(
                    numpy.float32
                ) / corrected_intensity_seg.astype(
                    numpy.float32
                )
        # If seg. vs itself, then just write a 1.0
        else:
            properties_row["ratio_{}_vs_seg".format(c)] = 1.0

    properties_row["info_fov"] = fov
    properties_row["info_frame"] = frame_number
    properties_row["info_seg_channel"] = enum
    properties_row["info_label"] = region.label
    properties_row["info_trench_number"] = z

    properties_row["geometry_Area"] = region.area
    properties_row["geometry_Orientation"] = region.orientation
    properties_row["geometry_Perimeter"] = region.perimeter

    # FIXME very small numbers --> math domain error.
    # Shouldn't the skimage library handle this issue and just return NaN?
    # The values might drop below zero due to floating-point errorm
    # and then a square root of a neg number causes the error.
    # Maybe we can wrap each of these in try clauses, and write NaN if they throw exceptions.
    properties_row["axis_length_major"] = region.major_axis_length
    properties_row["axis_length_minor"] = region.minor_axis_length

    properties_row["centroid_row"] = region.centroid[0]
    properties_row["centroid_col"] = region.centroid[1]

    properties_row["bounding_box_min_row"] = region.bbox[0]
    properties_row["bounding_box_min_col"] = region.bbox[1]
    properties_row["bounding_box_max_row"] = region.bbox[2]
    properties_row["bounding_box_max_col"] = region.bbox[3]

    ### Done iterating seg channels
    # Append the region information to the table
    properties_row.append()


# Input a list of pixel co-ordinates (as per skimage.measure.regionpprops),
# an associated image (2d numpy array) to read intensities from, and and a background signal value.
# Return the sum of the pixel intensities minus the product of the background value and the number of pixels.
def subtract_background_from_coords(coords, image, background):
    summed_intensity = image[coords[0], coords[1]].sum()
    total_background = len(coords) * background

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


# Merge the tables at the very end
def merge_tables(in_dir, identifier, parent_file, fov_dir, channels, seg_channels):
    Cell = make_cell_type(channels, seg_channels)

    h5file = tables.open_file(os.path.join(in_dir, parent_file), mode="r+")
    table_name = "measurements_thresholded_trench_{}".format(identifier)
    table_path = "/tables/{}".format(table_name)

    # Create a table for storing the properties of the entire trench
    # If table exists, delete it and start over
    if h5file.__contains__(table_path):
        node = h5file.get_node(table_path)
        node._f_remove()

    cell_properties_table_global = h5file.create_table(
        h5file.root.tables, table_name, Cell, "Measurements_thresholded_trench"
    )

    ### See: https://stackoverflow.com/questions/24891014/append-all-rows-from-one-table-to-another-using-pytables

    for filename in os.listdir(fov_dir):
        if filename.endswith("h5"):
            table_file = tables.open_file(os.path.join(fov_dir, filename), "r")
            table = table_file.get_node(table_path)

            for next_row in table.iterrows():
                cell_properties_table_global.append([next_row[:]])

            table.close()
            table_file.close()

    # Done copying trench properties
    cell_properties_table_global.flush()
    cell_properties_table_global.close()

    h5file.close()


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
        "-P", "--params-file", type=str, help="Segmentation parameters file"
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
    params_file_path = args.params_file

    # FIXME what's the right way to read these in?
    x_dimension = 30
    y_dimension = crop_bottom - crop_top

    ###

    # Load the ND2 metadata to get the channel names
    # NOTE do we need this?
    # channel_names = get_channel_names(in_dir, parent_file)

    params = read_params_file(params_file_path)

    # Let's just focus on fluor. for now.
    channels = parse_fluorescence_channel_names(params)

    seg_channels = parse_segmentation_channel_names(params)

    fov_dir = os.path.join(in_dir, "FOV")

    files = os.listdir(fov_dir)
    # DEBUG
    # files = ["FOV_0.h5"]

    ### Iterate all the files
    # NOTE: this could also be done outside of Python (e.g. in a Slurm / Bash script)
    # in an embarrasingly parallel fashion.

    print("Measuring cell properties...")
    start = time.time()

    # Run in parallel
    if num_cpu > 1:
        args = [
            (
                os.path.join(fov_dir, filename),
                identifier,
                channels,
                seg_channels,
                crop_top,
                crop_bottom,
                x_dimension,
                y_dimension,
                params,
            )
            for filename in files
        ]

        # TODO: weed out any files which do not end in h5

        with Pool(processes=num_cpu) as p:
            p.starmap(run_segmentation_analysis, args)

    # Do not run in parallel
    else:
        for filename in files:
            if filename.endswith("h5"):
                run_segmentation_analysis(
                    os.path.join(fov_dir, filename),
                    identifier,
                    channels,
                    seg_channels,
                    crop_top,
                    crop_bottom,
                    x_dimension,
                    y_dimension,
                    params,
                )

    ### Done looping all FOV
    end = time.time()
    print(end - start)
    print("Done measuring cell properties.")

    ### Merge the tables
    print("Merging tables...")
    start = time.time()

    # TODO: also store the BF/phase summed intensity information?
    merge_tables(in_dir, identifier, parent_file, fov_dir, channels, seg_channels)

    end = time.time()
    print(end - start)
    print("Done merging tables.")
