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


### Segment cells in trenches
###
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
                these_params[param_names[i + 1]] = e

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
def make_cell_type(channels):
    # Define the PyTables column types using a dictionary
    column_types = {
        "info_fov": tables.UInt16Col(),
        "info_frame": tables.UInt16Col(),
        "info_seg_channel": tables.EnumCol(channels, channels[0], "int8"),
        "info_label": tables.UInt16Col(),
        "info_trench_number": tables.UInt16Col(),
        "geometry_Area": tables.UInt16Col(),
        "geometry_Orientation": tables.Float32Col(),
        "geometry_Perimeter": tables.Float32Col(),
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
        column_types["total_intensity_{}".format(c)] = tables.UInt32Col()

    return type("Cell", (tables.IsDescription,), column_types)


# Return average pixel values in a block_size*block_size region
# Modified from https://github.com/scikit-image/scikit-image/blob/v0.16.1/skimage/filters/thresholding.py#L142
def calculate_local_mean(image, block_size):
    thresh_image = numpy.zeros(image.shape, "double")
    mask = 1.0 / block_size * numpy.ones((block_size,))
    # separation of filters to speedup convolution
    ndimage.convolve1d(image, mask, axis=0, output=thresh_image, mode="reflect")
    ndimage.convolve1d(thresh_image, mask, axis=1, output=thresh_image, mode="reflect")

    return thresh_image


# Numpy version of the segmentation algorithm
def run_segmentation_numpy(
    stack,
    stack_phase,
    threshold_fluor,
    threshold_phase,
    window_size,
    k,
    x_dimension,
    y_dimension,
    block_size,
):
    # Absolute thresholding
    # TODO: optionally apply local averaging to mitigate noisy pixels
    # stack_bimodal = (stack < threshold_fluor) * 0 + (stack >= threshold_fluor) * stack

    # Absolute thresholding after smoothing out noise using a local averaging of pixel intensities
    # in a block_sizexblock_size region.
    # Set pixels < threshold to 0, > threshold to their original value, but using the smoothed values for comparison.
    # This helps deal with lone pixels which are noisy & low intensity: prevent them from creating gaps in a cell.
    stack_local_mean = calculate_local_mean(stack, block_size)
    stack_bimodal = (stack_local_mean < threshold_fluor) * 0 + (
        stack_local_mean >= threshold_fluor
    ) * stack

    # Include phase information: discard bright pixels
    # NOTE trenches in phase have uneven brightness -- doesn't work so well
    stack_phase_bimodal = stack_phase < threshold_phase
    stack_bimodal = stack_bimodal * stack_phase_bimodal

    # Binary mask
    binary_bimodal = stack_bimodal > 0

    # Run an opening (erosion followed by dilation) to get rid of lone pixels and to smooth the cell contours
    # Mahotas requires looping over each z level
    for z in range(binary_bimodal.shape[2]):
        binary_bimodal[:, :, z] = mahotas.morph.open(binary_bimodal[:, :, z])
        #         binary_bimodal[:, :, z] = mahotas.morph.erode(binary_bimodal[:, :, z])
        binary_bimodal[:, :, z] = mahotas.morph.dilate(binary_bimodal[:, :, z])

    # And apply this operator's result back onto the bimodal_stack
    stack_bimodal = (binary_bimodal == 0) * 0 + (binary_bimodal == 1) * stack_bimodal

    # FIXME small object removal when??
    # Remove small objects
    for z in range(stack_bimodal.shape[2]):
        stack_bimodal = skimage.morphology.remove_small_objects(
            stack_bimodal, min_size=40, connectivity=2
        )

    # Apply Niblack method
    niblack = numpy.empty(stack.shape, dtype=numpy.int16)
    for z in range(stack.shape[2]):
        niblack[:, :, z] = skimage.filters.threshold_niblack(
            stack_bimodal[:, :, z], window_size=window_size, k=k
        )

    segmentation_mask = (stack > niblack) * binary_bimodal

    # FIXME small object removal when??
    # Remove small objects
    for z in range(segmentation_mask.shape[2]):
        stack_bimodal = skimage.morphology.remove_small_objects(
            segmentation_mask, min_size=100, connectivity=2
        )

    # Label the regions
    regions = numpy.empty(
        shape=(y_dimension, x_dimension, stack.shape[2]), dtype=numpy.int16
    )
    for z in range(stack.shape[2]):
        regions[:, :, z], n_regions = mahotas.label(segmentation_mask[:, :, z])

    # Multiply by -1 to make it a basin
    # And cast to signed integer because the skimage watershed algorithm requires it
    img = (stack * -1).astype(numpy.int16)

    return (regions, img, binary_bimodal)


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
    stack = (
        ((stack > background) * stack) + ((stack <= background) * background)
    ) - background

    return stack


# TODO: is there a way to merge some of this code with the whole-trench code? Reduce duplication...
# FIXME read ranges from the detected trenches table, extract from it the x_y dimensions & the stack size!!
def run_segmentation_analysis(
    fov_file_path,
    identifier,
    channel_names,
    crop_top,
    crop_bottom,
    x_dimension,
    y_dimension,
    params,
    fluorescence_channels,
    segmentation_channels,
    background_correction,
    bf_channel_name,
    minimum_cell_area,
    block_size,
):

    # For writing the segmentation masks
    filters = tables.Filters(complevel=1, complib="zlib")
    chunk_dimensions = (128, 128)

    # BF or Phase cutoff
    phase_cutoff = int(params["BF"]["threshold"])

    # Is there a way to pass in the function, if we can't pass in the type?
    # FIXME all channel names, or just the fluorescence channels?
    Cell = make_cell_type(channel_names)

    # Define the enumerable type, for storing in the table
    segmentation_enum = tables.Enum(segmentation_channels)

    # Open this FOV file
    fov_h5file = tables.open_file(fov_file_path, mode="r+")

    # Extract the FOV number from the file name
    base_name = os.path.basename(fov_file_path)
    (root, extension) = os.path.splitext(base_name)
    fov = int(root.split("_")[1])

    table_name = "measurements_segmented_cell_{}".format(identifier)

    # Create the properties table
    # Create properties tables for each FOV, then merge them together at the end
    if fov_h5file.__contains__("/tables/{}".format(table_name)):
        node = fov_h5file.get_node("/tables/{}".format(table_name))
        node._f_remove()

    properties_table = fov_h5file.create_table(
        "/tables", table_name, Cell, "Measurements_segmented_cell"
    )
    properties_row = properties_table.row

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

        ### Segmentation
        # For each segmentation channel,
        # segment & measure; write masks & measurements to table

        # Input channel name, return stack of trenches
        channel_to_stack = {}
        for c in channel_names:
            channel_to_stack[c] = load_stack_trenches(
                time_node._f_get_child(c),
                crop_top,
                crop_bottom,
                x_dimension,
                y_dimension,
                stack_size,
                ranges,
                background_correction,
            )

        # Segmentation on the whole stack of trenches
        for sc in segmentation_channels:
            (regions, basin, thresholded) = run_segmentation_numpy(
                channel_to_stack[sc],
                channel_to_stack[bf_channel_name],
                int(params[sc]["threshold"]),
                phase_cutoff,
                int(params[sc]["niblack_window_size"]),
                float(params[sc]["niblack_k"]),
                x_dimension,
                y_dimension,
                block_size,
            )

            # The enumerable type for table columns
            sc_to_enum = segmentation_enum[sc]

            # Watershed on each trench
            # Run measurements on regions & store info. in HDF5 table
            for n in range(regions.shape[2]):
                watershed_result = skimage.morphology.watershed(
                    image=basin[:, :, n],
                    markers=regions[:, :, n],
                    mask=thresholded[:, :, n],
                )

                # Store the segmentation result (the "mask") in the HDF5 file
                # NOTE even empty segmentation masks will be stored. Is this necessary?
                fov_h5file.create_carray(
                    "/masks/Frame_{}/Channel_{}".format(frame_number, sc),
                    "trench_{}".format(n),
                    obj=watershed_result,
                    title="{}/{} segmentation mask".format(sc, n),
                    chunkshape=chunk_dimensions,
                    filters=filters,
                    createparents=True,
                )

                # Were any cells actually detected?
                unique_ids = numpy.unique(watershed_result)

                if len(unique_ids) > 1:
                    properties = skimage.measure.regionprops(watershed_result)

                    for cell in properties:
                        if cell.area >= minimum_cell_area:
                            transposed = cell.coords.T

                            # Get the sum of all the pixel intensities of the cell region
                            for c in channel_names:
                                properties_row[
                                    "total_intensity_{}".format(c)
                                ] = channel_to_stack[c][:, :, n][
                                    transposed[0], transposed[1]
                                ].sum()

                            properties_row["info_fov"] = fov
                            properties_row["info_frame"] = frame_number
                            properties_row["info_seg_channel"] = sc_to_enum
                            properties_row["info_label"] = cell.label
                            properties_row["info_trench_number"] = n

                            properties_row["geometry_Area"] = cell.area
                            properties_row["geometry_Orientation"] = cell.orientation
                            properties_row["geometry_Perimeter"] = cell.perimeter

                            properties_row["axis_length_major"] = cell.major_axis_length
                            properties_row["axis_length_minor"] = cell.minor_axis_length

                            properties_row["centroid_row"] = cell.centroid[0]
                            properties_row["centroid_col"] = cell.centroid[1]

                            properties_row["bounding_box_min_row"] = cell.bbox[0]
                            properties_row["bounding_box_min_col"] = cell.bbox[1]
                            properties_row["bounding_box_max_row"] = cell.bbox[2]
                            properties_row["bounding_box_max_col"] = cell.bbox[3]

                            # Append the cell information to the table
                            properties_row.append()

    ### Done measuring trenches in this FOV, across all time points
    properties_table.flush()
    properties_table.close()

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
    stack = (
        ((stack > background) * stack) + ((stack <= background) * background)
    ) - background

    return stack


# Merge the tables at the very end
def merge_tables(in_dir, identifier, parent_file, fov_dir, channel_names):
    Cell = make_cell_type(channel_names)

    h5file = tables.open_file(os.path.join(in_dir, parent_file), mode="r+")
    table_name = "measurements_segmented_cell_{}".format(identifier)

    # Create a table for storing the properties of the entire trench
    # If table exists, delete it and start over
    if h5file.__contains__("/tables/{}".format(table_name)):
        node = h5file.get_node("/tables/{}".format(table_name))
        node._f_remove()

    cell_properties_table_global = h5file.create_table(
        h5file.root.tables, table_name, Cell, "Measurements_segmented_cell"
    )

    ### See: https://stackoverflow.com/questions/24891014/append-all-rows-from-one-table-to-another-using-pytables

    for filename in os.listdir(fov_dir):
        if filename.endswith("h5"):
            table_file = tables.open_file(os.path.join(fov_dir, filename), "r")
            table = table_file.get_node("/tables/{}".format(table_name))

            for next_row in table.iterrows():
                cell_properties_table_global.append([next_row[:]])

            table.close()
            table_file.close()

    # Done copying trench properties
    cell_properties_table_global.flush()
    cell_properties_table_global.close()


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
    background_correction = 100
    bf_channel_name = "BF"
    minimum_cell_area = 50
    block_size = 9

    ###

    # Load the ND2 metadata to get the channel names
    # NOTE do we need this?
    channel_names = get_channel_names(in_dir, parent_file)

    params = read_params_file(params_file_path)

    fluorescence_channels = parse_fluorescence_channel_names(params)

    segmentation_channels = parse_segmentation_channel_names(params)

    fov_dir = os.path.join(in_dir, "FOV")

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
                channel_names,
                crop_top,
                crop_bottom,
                x_dimension,
                y_dimension,
                params,
                fluorescence_channels,
                segmentation_channels,
                background_correction,
                bf_channel_name,
                minimum_cell_area,
                block_size,
            )
            for filename in os.listdir(fov_dir)
        ]

        # TODO: weed out any files which do not end in h5

        with Pool(processes=num_cpu) as p:
            p.starmap(run_segmentation_analysis, args)

    # Do not run in parallel
    else:
        for filename in os.listdir(fov_dir):
            if filename.endswith("h5"):
                run_segmentation_analysis(
                    os.path.join(fov_dir, filename),
                    identifier,
                    channel_names,
                    crop_top,
                    crop_bottom,
                    x_dimension,
                    y_dimension,
                    params,
                    fluorescence_channels,
                    segmentation_channels,
                    background_correction,
                    bf_channel_name,
                    minimum_cell_area,
                    block_size,
                )

    ### Done looping all FOV
    end = time.time()
    print(end - start)
    print("Done measuring cell properties.")

    ### Merge the tables
    print("Merging tables...")
    start = time.time()

    # TODO: also store the BF/phase summed intensity information?
    merge_tables(in_dir, identifier, parent_file, fov_dir, channel_names)

    end = time.time()
    print(end - start)
    print("Done merging tables.")
