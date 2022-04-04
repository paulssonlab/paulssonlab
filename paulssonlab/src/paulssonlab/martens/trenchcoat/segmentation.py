#!/usr/bin/env python3

import skimage
import numpy
import os
import pathlib
import tables
from multiprocessing import Pool
from tqdm import tqdm
import pandas

import algorithms
import algo_sharp

# import arrayfire_algorithms
from properties import (
    # write_properties_to_table,
    # merge_tables,
    make_cell_type,
    # subtract_background_from_coords,
    init_properties_dict,
    write_properties_to_table_from_df,
    get_max_length,
    add_properties,
)
from metadata import get_metadata
from params import read_params_file
from napari_browse_hdf5 import parse_corrections_settings, sub_bg_no_underflow

"""
Perform cell segmentation & measure fluorescence intensities in microscope images, with support for sub-regions (e.g. "trenches").

TODO:
1. Pass in min. region size. Use pixel_microns to help convert from true area in microns^2 to region size in pixels.
   Makes this parameter magnification-independent.

2. Are the masks being symbolically linked correctly in the HDF5 files?
"""


def run_segmentation_analysis_regions(
    in_file,
    name,
    fov,
    frame,
    z_level,
    seg_params,
    channels,
    out_dir_masks,
    algo_dict,
    regions_file,
    max_len_filenames,
    max_len_seg_channels,
    corrections_dict,
):
    """
    Read region co-ordinates from an HDF5 file (with support for multiple sets of regions per image).
    Then call the code to do the actual segmentation work.
    """
    # Store results here
    # For each entry in the list, a dict -> numpy arrays
    results = []

    # H5file with image data
    h5file = tables.open_file(in_file, mode="r")
    z_node = h5file.get_node(
        "/Images/{}/FOV_{}/Frame_{}/Z_{}".format(name, fov, frame, z_level)
    )

    # H5file with a table of trench co-ordinates
    h5file_reg = tables.open_file(regions_file, mode="r")

    # Read in the array of regions from a table
    # 1. Query the table for the relevant rows
    trench_table = h5file_reg.get_node("/trench_coords")
    query = """(info_file == name) & (info_fov == fov) & (info_frame == frame) & (info_z_level == z_level)"""
    this_node = trench_table.read_where(query)

    # 2. Trench rows: unique entries in the "info_row_number" column
    df = pandas.DataFrame(this_node)
    rows = df["info_row_number"].unique()

    # 3. Compile the rows into a numpy array & run segmentation
    for row_number in rows:
        # a. Grab all the relevant trenches
        trenches = df[(df["info_row_number"] == row_number)]
        # FIXME is it imperative to sort the rows in ascending order?
        # (they must be ascending; but maybe they will always be sorted, since that's the order they were written?)

        # b. Convert this set of trenches into a single regions array
        num_trenches = len(trenches)
        regions = numpy.empty((num_trenches, 4), dtype=numpy.uint16)
        for i, r in enumerate(trenches.itertuples()):
            regions[i] = [r.min_row, r.min_col, r.max_row, r.max_col]

        # c. Use the array of bounding boxes to extract regions from the image & create an image stack
        (stack, ch_to_index) = make_ch_to_img_stack_regions(
            h5file, z_node, channels, regions
        )

        # Apply flat-fielding corrections?
        if corrections_dict is not None:
            stack = run_flatfielding(stack, ch_to_index, regions, corrections_dict)

        # d. Run the analysis on the given image stack
        # NOTE would it make sense to also pass in the region info,
        # so that the cell centroids can also be recorded with respect
        # to the entire image, and not just the trench bounding box?
        # (when there are no regions, would have to be left empty,
        # (or identical)).

        # Contains results, 1 per seg channel
        results_this_row = run_segmentation_analysis(
            in_file,
            name,
            fov,
            frame,
            z_level,
            seg_params,
            channels,
            out_dir_masks,
            algo_dict,
            row_number,
            stack,
            ch_to_index,
            True,
            max_len_filenames,
            max_len_seg_channels,
            corrections_dict,
        )

        # Contatenate the lists from each row of trenches
        # NOTE alternatively, we could pass in the list down a level,
        # and append to the same one always.
        for entry in results_this_row:
            results.append(entry)

    h5file_reg.close()
    h5file.close()

    return results


def run_flatfielding(stack, ch_to_index, regions, corrections_dict):
    """
    Given a stack + channel-to-index,
    and given a set of corrections,
    apply the corrections to the data.

    Corrections include:
    1. Subtracting "dark" image (camera noise), or a fixed value (e.g. 100), taking care not to underflow
    2. Dividing the sample image by the flatfielding image (i.e. fluorescent dyes)

    Corrections is a dict which maps between channels & correction matrices. Same convention as used in
    napari_browse_hdf5.py.

    TODO registration? (beads correction)

    TODO what if there are regions? Then need to slice out each specific regional
    area from the correction images, then do the corrections, one region at a time.
    """
    if corrections_dict["camera_noise"] is not None:
        camera_noise = corrections_dict["camera_noise"]
    else:
        # A default value of 100.0
        camera_noise = 100.0

    # No Regions ~ could also check if regions is None
    if stack.shape[2] == 1:
        # 1. Camera Noise
        for (channel, index) in ch_to_index.items():
            stack[..., 0, index] = sub_bg_no_underflow(
                stack[..., 0, index], camera_noise
            )

        # 2. Flatfielding
        # Optional for each channel
        # If there is nothing, then do nothing!
        if corrections_dict["flatfield"] is not None:
            for (channel, index) in ch_to_index.items():
                if corrections_dict["flatfield"][channel] is not None:
                    stack[..., 0, index] = numpy.divide(
                        stack[..., 0, index], corrections_dict["flatfield"][channel]
                    )

    # Regions
    else:
        if corrections_dict["flatfield"] is not None:
            for (channel, index) in ch_to_index.items():
                if corrections_dict["flatfield"][channel] is not None:
                    # min_row, min_col, max_row, max_col
                    for j, r in enumerate(regions):
                        stack[..., j, index] = image_node[
                            r[0] : r[2], r[1] : r[3]
                        ].astype(numpy.float32, casting="safe", copy=False)
                        # Add warning if values become negative? How to tell without checking all pixels? (could be expensive)
                        # could check if min val is > 0 after running operation.
                        stack[..., j, index] -= camera_bias[channel]
                        stack[..., j, index] /= f_field[channe][
                            r[0] : r[2], r[1] : r[3]
                        ]

    return stack


def run_segmentation_analysis_no_regions(
    in_file,
    name,
    fov,
    frame,
    z_level,
    seg_params,
    channels,
    out_dir_masks,
    algo_dict,
    max_len_filenames,
    max_len_seg_channels,
    corrections_dict,
):
    """
    Open h5file with images, load without regions. Call segmentation.
    """
    # Node in the input file, which may contain multiple channel images
    h5file = tables.open_file(in_file, mode="r")
    z_node = h5file.get_node(
        "/Images/{}/FOV_{}/Frame_{}/Z_{}".format(name, fov, frame, z_level)
    )

    (stack, ch_to_index) = make_ch_to_img_stack_no_regions(h5file, z_node, channels)
    h5file.close()

    # Apply flat-fielding corrections?
    if corrections_dict is not None:
        stack = run_flatfielding(stack, ch_to_index, None, corrections_dict)

    # No regions, no regions_set_number -> set row to 0
    results = run_segmentation_analysis(
        in_file,
        name,
        fov,
        frame,
        z_level,
        seg_params,
        channels,
        out_dir_masks,
        algo_dict,
        0,
        stack,
        ch_to_index,
        False,
        max_len_filenames,
        max_len_seg_channels,
        corrections_dict,
    )

    # No regions, no rows, so pass along the stack, ch_to_index,
    # and the global segmentation results for each seg channel:
    # the masks, & associated data (name, fov, etc.)
    return results


def run_segmentation_analysis(
    in_file,
    name,
    fov,
    frame,
    z_level,
    seg_params,
    channels,
    out_dir_masks,
    algo_dict,
    row_number,
    stack,
    ch_to_index,
    has_regions,
    max_len_filenames,
    max_len_seg_channels,
    corrections_dict,
):
    """
    Create new HDF5 files for writing masks, measurements table Load all
    channel images within a given File / FOV / Frame / Z-level, and extract
    regions using region co-ordinate Call function to do segmentation,
    measurements, and writing results to disk.

    This function is shared by both the regions and no_regions codepaths (where they converge).
    """

    # Create directory structure to store the masks
    # NOTE Must set exist_ok=True to allow multiple frames or regions,
    # but now we have the possibility that someone running the same analysis multiple times
    # will just keep on appending data over & over!
    # Thus it is essential to manually delete the directory structure if re-running an analysis.
    pathlib.Path("{}/{}/FOV_{}/".format(out_dir_masks, name, fov)).mkdir(
        parents=True, exist_ok=True
    )
    h5file_masks = tables.open_file(
        "{}/{}/FOV_{}/Frame_{}.h5".format(out_dir_masks, name, fov, frame), mode="a"
    )

    # This function calls the segmentation code & cell measurements code, and writes the masks to disk.
    # It also returns a list of batched results, 1 entry for each seg channel
    # Each batch is a dict of numpy arrays, one for each column of data.
    results = write_masks(
        h5file_masks,
        name,
        fov,
        frame,
        z_level,
        seg_params,
        stack,
        ch_to_index,
        algo_dict,
        row_number,
        max_len_filenames,
        max_len_seg_channels,
        channels,
    )

    # Done!
    h5file_masks.close()

    return results


def write_masks(
    h5file_masks,
    name,
    fov,
    frame,
    z_level,
    seg_params,
    stack,
    ch_to_index,
    algo_dict,
    row_number,
    max_len_filenames,
    max_len_seg_channels,
    channels,
):
    """
    Compute masks using specified segmentation algorithm.

    Write masks & measurements to HDF5 files. Analyze regions within
    images. (e.g. trenches, cropped images...)

    Return a list of dict of numpy arrays.
    seg channel -> dict
    dict -> arrays of cell measurements
    """
    batches = []

    for sc in seg_params.keys():
        # Retrieve the segmentation function
        function = algo_dict[sc]

        # Calculate the mask(s)
        masks = function(stack, ch_to_index, seg_params[sc])

        # Now, we know exactly how many masks there are total,
        # across all regions. This allows us to set a length
        # for this batch of entries.
        # FIXME the arrays might be larger than needed,
        # which causes pandas to convert to float to avoid empty values.
        # Is this because we mis-count the number of cells?
        # A possible fix includes converting back from python lists,
        # but this is a pain because of the dtypes.
        # It *should* be possible to exactly predict the sizes!!
        total_num_cells = 0

        # If we determine the number of unique pixel labels in each region,
        # then we can pre-calculate the total number of cells.
        for region_number in range(masks.shape[2]):
            mask = masks[..., region_number]

            # Is it _not_ all zero?
            if mask.any():
                len_uni = len(numpy.unique(mask))

                # Minus 1, because zero (background) are also values
                # But only make this correction if there are any background
                # values.
                if mask.min() == 0:
                    len_uni -= 1

                total_num_cells += len_uni

        # Initialize the structure
        if masks.shape[2] > 1:
            has_regions = True
        else:
            has_regions = False

        this_batch = init_properties_dict(
            channels,
            has_regions,
            total_num_cells,
            max_len_filenames,
            max_len_seg_channels,
        )

        # Write the masks
        global_index = 0
        for region_number in range(masks.shape[2]):
            mask = masks[..., region_number]
            # Write to disk
            # TODO delay writing to disk alongside the table?
            # Change how we do this, make a single large H5file instead
            # of breaking into chunks?
            h5file_masks.create_carray(
                "/Z_{}/{}/{}".format(z_level, sc, row_number),
                "region_{}".format(region_number),
                obj=mask,
                createparents=True,
                chunkshape=(128, 128),
                filters=tables.Filters(complevel=1, complib="zlib"),
            )

            # Compute the properties
            # NOTE / TODO: future versions of skimage will allow spitting out a properties object all at once,
            # rather than lazily calculating them one at a time.
            # TODO skip these steps if we can determine the mask is empty?
            cell_properties = skimage.measure.regionprops(mask)

            for cell in cell_properties:
                add_properties(
                    this_batch,
                    mask,
                    name,
                    fov,
                    frame,
                    z_level,
                    row_number,
                    region_number,
                    sc,
                    stack,
                    ch_to_index,
                    cell,
                    global_index,
                )

                global_index += 1

        batches.append(this_batch)

    # A list of batched results, 1 entry for each seg channel
    # Each batch is a dict of numpy arrays, one for each column of data.
    return batches


def make_ch_to_img_stack_regions(h5file, z_node, channels, regions):
    """
    Input an H5file, Z-node, channel names, region locations
    Returns a 4-D stack of all channels & regions
    0. channel 1. regions 2. X 3. Y

    Ranges is a numpy array with (min_row, min_col, max_row, max_col) for each region.
    NOTE: Pass in the node reference, not the whole image, into this sub-routine.
    Then, when using the slice notation below, will take advantage of
    the chunking so as to not load the entire image.

    Forcibly load all images as single-precision floats.

    NOTE: the skimage libraries technically follow the convention that images are between 0 and 1.
    This requires resampling data to that range, and then converting it back to the original range.
    Need to check for which algorithms this matters (Niblack? Otsu?).

    TODO: pass in the image loading function & image variable,
    to allow different kinds of image loading (not just HDF5).
    e.g. For HDF5, we would pass in the image_node.
    Might not even need to pass in a function, if all images
    support bracketed indexing.
    How to handle multiple images though?
    """

    # Assume all regions have the same dimensions.
    # Use the zeroth element to calculate the dimensions.
    x_dimension = regions[0, 2] - regions[0, 0]
    y_dimension = regions[0, 3] - regions[0, 1]

    ch_to_index = {}
    # Initialize the array as empty, with the 'F' ordering, and then to load it x/y/r/c
    stack = numpy.empty(
        (x_dimension, y_dimension, regions.shape[0], len(channels)),
        dtype=numpy.float32,
        order="F",
    )

    for i, ch in enumerate(channels):
        ch_to_index[ch] = i

        image_node = h5file.get_node(z_node, ch)

        # TODO: here, add optional camera bias & flat-field corrections to the image, before
        # loading the regions.
        # Biases: 1 per channel
        # Flat-field corrections: 1 per channel (could re-use the same ones, if necessary)
        # Allow entries to be none? (and skippable?)
        # OR do this *after* loading, to take advantage of chunking.
        #
        # # min_row, min_col, max_row, max_col
        # for j, r in enumerate(regions):
        #     stack[..., j, i] = image_node[r[0] : r[2], r[1] : r[3]].astype(numpy.float32, casting="safe", copy=False)
        #     # Add warning if values become negative? How to tell without checking all pixels? (could be expensive)
        #     # could check if min val is > 0 after running operation.
        #     stack[..., j, i] -= camera_bias[channel]
        #     stack[..., j, i] /= f_field[channe][r[0] : r[2], r[1] : r[3]]

        # min_row, min_col, max_row, max_col
        for j, r in enumerate(regions):
            # NOTE Does copy flag actually make a difference?
            stack[..., j, i] = image_node[r[0] : r[2], r[1] : r[3]].astype(
                numpy.float32, casting="safe", copy=False
            )

    return (stack, ch_to_index)


def make_ch_to_img_stack_no_regions(h5file, z_node, channels):
    """
    Input an H5file, Z-node, channel names,
    Returns a 4-D stack of all channels & a dummy region of dimension 1
    (for compatibility with the regions version)
    0. channel 1. dummy region 2. X 3. Y

    Forcibly load all images as single-precision floats.

    NOTE: Pass in the node reference, not the whole image, into this sub-routine.
    Then, when using the slice notation below, will take advantage of
    the chunking so as to not load the entire image.

    NOTE: the skimage libraries technically follow the convention that images are between 0 and 1.
    This requires resampling data to that range, and then converting it back to the original range.
    Need to check for which algorithms this matters (Niblack? Otsu?).
    """

    ch_to_index = {}

    # Have to process the zeroth image to get its dimensions, before looping over the remaining ones
    zeroth_channel = channels[0]
    ch_to_index[zeroth_channel] = 0

    # NOTE Does copy flag actually make a difference?
    zeroth_img = (
        h5file.get_node(z_node, zeroth_channel)
        .read()
        .astype(numpy.float32, casting="safe", copy=False)
    )

    # Initalize the empty ndarray
    stack = numpy.empty(
        (zeroth_img.shape[0], zeroth_img.shape[1], 1, len(channels)),
        dtype=numpy.float32,
        order="F",
    )

    # Pad an extra axis
    zeroth_img = zeroth_img[..., numpy.newaxis]
    stack[..., 0] = zeroth_img

    # Loop over the remaining images (channels)
    iter_channels = enumerate(channels)
    # Skip the zeroth one, because we already handled it
    next(iter_channels)
    for i, ch in iter_channels:
        ch_to_index[ch] = i

        # NOTE Does copy flag actually make a difference?
        image = (
            h5file.get_node(z_node, ch)
            .read()
            .astype(numpy.float32, casting="safe", copy=False)
        )

        image = image[..., numpy.newaxis]
        stack[..., i] = image

    return (stack, ch_to_index)


def link_files(in_dir, file_name):
    """Link the masks or tables files into a single HDF5 file."""

    out_file = os.path.join(in_dir, "{}.h5".format(file_name))
    h5file_top = tables.open_file(out_file, mode="w")
    # ND2 file directories
    for nd2_entry in os.scandir(in_dir):
        if nd2_entry.is_dir():
            # H5 file for each ND2 file
            nd2_h5_outfile = os.path.join(in_dir, "{}.h5".format(nd2_entry.name))
            h5file_nd2file = tables.open_file(
                nd2_h5_outfile, mode="w", title=nd2_entry.name
            )

            # FOV directories
            for fov_entry in os.scandir(nd2_entry):
                if fov_entry.is_dir():
                    h5file_fov = tables.open_file(
                        "{}.h5".format(fov_entry.path), mode="w", title=fov_entry.name
                    )

                    # The individual H5 files, for each time frame (containing Z/channel stacks)
                    for frame_entry in os.scandir(fov_entry):
                        if frame_entry.is_file():
                            (name, extension) = os.path.splitext(
                                os.path.basename(frame_entry.name)
                            )
                            if extension == ".h5":
                                h5file_fov.create_external_link(
                                    "/",
                                    "{}".format(name),
                                    "{}/{}:/".format(fov_entry.name, frame_entry.name),
                                )

                    # Done with this FOV
                    h5file_fov.close()

                    # And link back to the h5_nd2 file
                    h5file_nd2file.create_external_link(
                        "/",
                        "{}".format(fov_entry.name),
                        "{}/{}.h5:/".format(nd2_entry.name, fov_entry.name),
                    )

            # Done with this file
            h5file_nd2file.close()

            # And link back to the top file
            h5file_top.create_external_link(
                "/{}".format(file_name),
                "{}".format(nd2_entry.name),
                "{}.h5:/".format(nd2_entry.name),
                createparents=True,
            )

    h5file_top.close()


def main_segmentation_function(out_dir, in_file, num_cpu, params_file, regions_file):
    """
    Input the HDF5 file, an out directory, how many CPUs, a YAML file with
    parameters, and an HDF5 file with region (e.g. trench) co-ordinates.

    Run cell segmentation & write masks and measurements to HDF5.
    """
    # Dir containing new HDF5 files, write results to
    out_dir_masks = os.path.join(out_dir, "MASKS")
    out_dir_tables = os.path.join(out_dir, "TABLES")

    # Create the masks & tables directories
    pathlib.Path(out_dir_masks).mkdir(parents=True, exist_ok=True)
    pathlib.Path(out_dir_tables).mkdir(parents=True, exist_ok=True)

    # Segmentation parameters, in YAML
    # and optional flatfielding & camera noise corrections parameters
    # TODO: verify that the channels specified in the params match the available channels in the files?
    # TODO: if an algorithm doesn't require any parameters, what do we do?
    params = read_params_file(params_file)

    # HDF5 file with images & metadata
    in_file = os.path.join(in_file, "data.h5")
    h5file = tables.open_file(in_file, mode="r")

    # Loop the nodes immediately under Images to get the file names
    file_names = [i._v_name for i in h5file.list_nodes("/Images")]
    # Max len is useful for setting column names with fixed string sizes
    max_len_filenames = get_max_length(file_names)

    # Get channels from one of the files
    # Assume that all files have the same channels: otherwise, cannot process them simultaneously!
    node = h5file.get_node("/Metadata/{}".format(file_names[0]))()
    channels = get_metadata(node)["channels"]
    channels = [c.decode("utf-8") for c in channels]

    # Iterate the files & images
    with Pool(processes=num_cpu) as pool:
        # Make a dict of algorithms
        algo_to_func = {
            "whole_trench": algorithms.measure_whole_trench,
            "threshold": algorithms.run_single_threshold,
            "niblack": algorithms.run_niblack_segmentation,
            "fluor_phase": algorithms.run_fluor_phase_segmentation,
            "fluor_sharpen": algo_sharp.run_fluor_sharp_segmentation
            #'dual_threshold'    : algorithms.run_dual_thresholding,
            #'niblack_phase_gpu' : algorithms.run_segmentation_GPU
        }

        # Input a channel, output an algorithm
        algo_dict = {}
        for ch, p in params["segmentation"].items():
            algo_dict[ch] = algo_to_func[p["algorithm"]]

        # Max len is useful for setting column names with fixed string sizes
        max_len_seg_channels = get_max_length(algo_dict.keys())

        # Loop the HDF5 file to count the num. of image stacks to process,
        # and use this value to initialize the progress bar.
        total = 0
        for n in h5file.iter_nodes(h5file.root.Images):
            metadata_node = h5file.get_node("/Metadata/{}".format(n._v_name))()
            metadata = get_metadata(metadata_node)
            total += (
                len(metadata["fields_of_view"])
                * len(metadata["frames"])
                * len(metadata["z_levels"])
            )

        # Parse the optional flatfielding, camera noise parameters
        try:
            corr_p = params["corrections"]
            corrections_dict = parse_corrections_settings(corr_p)
        except:
            print("No flatfielding params. detected.")
            corrections_dict = None

        # Analyze the images & write the masks
        print("Computing masks & measuring properties...")

        pbar = tqdm(total=total, desc="Channel stack")

        # Each task sent to the pool returns results back here
        # Currently, masks are written to disk by the processes themselves,
        # but cell measurements are passed back (in memory) to here.
        # The parent process therefore has to write these to the HDF5 table,
        # which is shared across all FOVs etc.
        # The idea is that the callback function prevents simultaneous writes
        # to the same HDF5 table, acting as a synchronizer between processes.
        # NOTE unclear if it's required to define within here,
        # so that it can "inherit" its parental scope? (Kind of ugly...)
        def write_tables(results):
            # Each list represents a column for a particular property
            # Each list should have exactly the same length (1 entry per cell).
            # Iterate the lists in parallel, & write them to the HDF5 table, one at a time.
            # OR convert this whole thing to a dataframe from_dict, and then write to HDF5.
            # Problem with from_dict: how to enforce dtype on a per-column basis?
            # https://stackoverflow.com/questions/42165705/assign-dtype-with-from-dict
            # https://github.com/pandas-dev/pandas/issues/14655
            # OR we first convert each list in the dict to a numpy array,
            # then assume that Pandas will know what to do?
            # see: dataframe_conversion.py

            # Have a separate batch of results for each row of trenches,
            # and for each segmentation channel (e.g.: 4 batches if
            # 2 rows & 2 seg channels).

            # Convert each batch from a dict of numpy arrays to a dataframe
            list_of_df = []
            for r in results:
                result_to_dict = pandas.DataFrame.from_dict(r)
                list_of_df.append(result_to_dict)

            # Concatenate them into a single dataframe
            merged_results = pandas.concat(list_of_df, axis=0)

            # Then, iterate the rows. Using a dataframe makes it syntactically less tedious
            # to grab a "row" of data and pass it into a subroutine.
            # Otherwise, we must grab every single possible column and pass it in,
            # which is weird because the number of columns isn't static (depends on
            # the number of seg channels, etc.), or we must copy the row into yet another dict.
            for cell_number in range(len(merged_results)):
                cell = merged_results.iloc[cell_number]
                write_properties_to_table_from_df(
                    cell, table_measurements.row, merged_results.columns
                )

            # NOTE it might also be possible to delay writing the masks until here.
            # However, this would require bubbling them up alongside the measurements.
            # Major advantage would be if we want to write masks to shared files.

            # Finally, pdate the progress bar
            pbar.update()

        # For printing error messages originating from sub-processes
        def error_callback(*e):
            print(e)

        # Create the table for cell measurments. Shared for all FOVs, frames etc.
        h5file_tables = tables.open_file(
            "{}/tables.h5".format(out_dir_tables), mode="w"
        )

        # Regions
        if regions_file:
            # allow extra columns for mother machine data.
            Cell = make_cell_type(channels, params.keys(), file_names, True)
            table_measurements = h5file_tables.create_table(
                "/",
                "measurements",
                Cell,
                createparents=True,
                filters=tables.Filters(complevel=1, complib="zlib"),
            )

            # Loop again, & perform the analysis
            for n in h5file.iter_nodes(h5file.root.Images):
                metadata_node = h5file.get_node("/Metadata/{}".format(n._v_name))()
                metadata = get_metadata(metadata_node)

                # TODO this would be a good place to calculate conversion from area
                # to pixels, metadata['using pixel_microns'], for remove_small_objects.

                for fov in metadata["fields_of_view"]:
                    for frame in metadata["frames"]:
                        for z_level in metadata["z_levels"]:
                            func_args = [
                                in_file,
                                n._v_name,
                                fov,
                                frame,
                                z_level,
                                params["segmentation"],
                                channels,
                                out_dir_masks,
                                algo_dict,
                                regions_file,
                                max_len_filenames,
                                max_len_seg_channels,
                                corrections_dict,
                            ]
                            pool.apply_async(
                                run_segmentation_analysis_regions,
                                func_args,
                                callback=write_tables,
                                error_callback=error_callback,
                            )
                            ##DEBUG
                            # results = run_segmentation_analysis_regions(*func_args)
                            # write_tables(results)

        # No regions
        else:
            # Don't bother with the extra columns for mother machine data (the width_intensity & Max_Width_Area).
            Cell = make_cell_type(channels, params.keys(), file_names, False)
            table_measurements = h5file_tables.create_table(
                "/",
                "measurements",
                Cell,
                createparents=True,
                filters=tables.Filters(complevel=1, complib="zlib"),
            )

            # Loop again, & perform the analysis
            for n in h5file.iter_nodes(h5file.root.Images):
                metadata_node = h5file.get_node("/Metadata/{}".format(n._v_name))()
                metadata = get_metadata(metadata_node)

                # TODO this would be a good place to calculate conversion from area
                # to pixels, metadata['using pixel_microns'], for remove_small_objects.

                for fov in metadata["fields_of_view"]:
                    for frame in metadata["frames"]:
                        for z_level in metadata["z_levels"]:
                            func_args = [
                                in_file,
                                n._v_name,
                                fov,
                                frame,
                                z_level,
                                params["segmentation"],
                                channels,
                                out_dir_masks,
                                algo_dict,
                                max_len_filenames,
                                max_len_seg_channels,
                                corrections_dict,
                            ]

                            pool.apply_async(
                                run_segmentation_analysis_no_regions,
                                func_args,
                                callback=write_tables,
                                error_callback=error_callback,
                            )
                            ##DEBUG
                            # results = run_segmentation_analysis_no_regions(*func_args)
                            # write_tables(results)

        pool.close()
        pool.join()

        pbar.close()
        table_measurements.flush()
        h5file_tables.close()

    h5file.close()

    print("Done computing masks & measuring properties.")

    # Link all the individual masks & properties files into respective H5 file,
    # to make it easier to iterate them.
    print("Linking masks...")
    link_files(out_dir_masks, "masks")

    # Write the segmentation params YAML to the masks h5file, so that the params can be referenced later on
    print("Saving YAML parameters...")
    params_file_handle = open(params_file)
    params_text = params_file_handle.read()
    params_file_handle.close()

    h5file_masks = tables.open_file(os.path.join(out_dir_masks, "masks.h5"), "r+")
    h5file_masks.create_array(
        "/Parameters",
        "seg_params.yaml",
        obj=numpy.array(params_text),
        createparents=True,
    )
    h5file_masks.close()
