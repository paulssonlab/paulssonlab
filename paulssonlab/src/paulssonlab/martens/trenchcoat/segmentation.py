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

# import arrayfire_algorithms
from properties import (
    write_properties_to_table,
    merge_tables,
    make_cell_type,
    subtract_background_from_coords,
)
from metadata import get_metadata
from params import read_params_file

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
    out_dir_tables,
    algo_dict,
    file_names,
    regions_file,
):
    """
    Read region co-ordinates from an HDF5 file (with support for multiple sets of regions per image).
    Then call the code to do the actual segmentation work.
    """
    # H5file with image data
    h5file = tables.open_file(in_file, mode="r")
    z_node = h5file.get_node(
        "/Images/{}/FOV_{}/Frame_{}/Z_{}".format(name, fov, frame, z_level)
    )

    # H5file with arrays denoting region co-ordinates
    # There can be multiple such arrays, e.g. if there are 2 rows of trenches,
    # then each row gets its own array, and each array specifies individual trenches 1 at a time.
    # NOTE now, no need to look at the arrays anymore, because all the info for
    # the entire dataset are compiled into a queryable table in advance.
    h5file_reg = tables.open_file(regions_file, mode="r")

    # Read in the array of regions from a table
    # 1. Query the table for the relevant rows
    trench_table = h5file_reg.get_node("/trench_coords")
    this_node = trench_table.read_where("""(info_file == name) & (info_fov == fov)""")
    # FIXME Don't query frame & z_level, since for now we are just using the first ones across all! ??
    # (info_frame == frame) & (info_z_level == z_level)

    # 2. Trench rows: unique entries in the "info_row_number" column
    df = pandas.DataFrame(this_node)
    region_sets = df["info_row_number"].unique()

    # 3. Compile the rows into a numpy array & run segmentation
    for region_set_number in region_sets:
        # a. Grab all the relevant trenches
        trenches = df[df["info_row_number"] == region_set_number]
        # FIXME is it imperative to sort the rows in ascending order?
        # (they must be ascending; but maybe they will always be sorted, since that's the order they were written?)

        # b. Convert this set of trenches into a single regions array
        # NOTE This was a first attempt using a dim=4 column, but that doesn't work with pandas.
        # Instead, have to merge 4 separate columns: min_row, min_col, max_row, max_col
        # regions = [i for i in trenches["bounding_box"]]
        # regions = numpy.array(regions)

        # Have a subset of dataframe rows which belong to trenches in this fov etc. & row of trenches
        # Take the 4 bbox coords for each trench & convert into a 2d array:
        # 1st dim is trench number, second dim is 4 bbox coords.
        regions = []
        for r in trenches.itertuples():
            regions.append([r.min_row, r.min_col, r.max_row, r.max_col])

        regions = numpy.array(regions)

        # c. Use the array of bounding boxes to extract regions from the image & create an image stack
        (stack, ch_to_index) = make_ch_to_img_stack_regions(
            h5file, z_node, channels, regions
        )

        # d. Run the analysis on the given image stack
        run_segmentation_analysis(
            in_file,
            name,
            fov,
            frame,
            z_level,
            seg_params,
            channels,
            out_dir_masks,
            out_dir_tables,
            algo_dict,
            file_names,
            region_set_number,
            stack,
            ch_to_index,
        )

    h5file_reg.close()
    h5file.close()


def run_segmentation_analysis_no_regions(
    in_file,
    name,
    fov,
    frame,
    z_level,
    seg_params,
    channels,
    out_dir_masks,
    out_dir_tables,
    algo_dict,
    file_names,
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

    # No regions, no regions_set_number -> set to 0
    run_segmentation_analysis(
        in_file,
        name,
        fov,
        frame,
        z_level,
        seg_params,
        channels,
        out_dir_masks,
        out_dir_tables,
        algo_dict,
        file_names,
        0,
        stack,
        ch_to_index,
    )


def run_segmentation_analysis(
    in_file,
    name,
    fov,
    frame,
    z_level,
    seg_params,
    channels,
    out_dir_masks,
    out_dir_tables,
    algo_dict,
    file_names,
    region_set_number,
    stack,
    ch_to_index,
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

    # Create directory structure to store the tables
    pathlib.Path("{}/{}/FOV_{}/".format(out_dir_tables, name, fov)).mkdir(
        parents=True, exist_ok=True
    )
    h5file_tables = tables.open_file(
        "{}/{}/FOV_{}/Frame_{}.h5".format(out_dir_tables, name, fov, frame), mode="a"
    )
    Cell = make_cell_type(channels, seg_params.keys(), file_names)

    # NOTE If multiple regions, don't try creating pre-existing table:
    # fix this by wrapping table creation in a try statement, and if that fails,
    # then open the existing table.
    try:
        table = h5file_tables.create_table(
            "/", "measurements", Cell, createparents=True
        )
    except:
        table = h5file_tables.get_node("/measurements")

    # This function calls the segmentation code & cell measurements code, and writes the results to disk.
    write_masks_tables(
        h5file_masks,
        h5file_tables,
        name,
        fov,
        frame,
        z_level,
        table.row,
        seg_params,
        stack,
        ch_to_index,
        algo_dict,
        region_set_number,
    )

    # Done!
    table.flush()
    h5file_masks.close()
    h5file_tables.close()


def write_masks_tables(
    h5file_masks,
    h5file_tables,
    name,
    fov,
    frame,
    z_level,
    row,
    seg_params,
    stack,
    ch_to_img,
    algo_dict,
    region_set_number,
):
    """
    Compute masks using specified segmentation algorithm.

    Write masks & measurements to HDF5 files. Analyze regions within
    images. (e.g. trenches, cropped images...)
    """

    for sc in seg_params.keys():
        # Calculate the mask(s)
        masks = algo_dict[sc](stack, ch_to_img, seg_params[sc])

        # Write the masks & tables
        for region_number in range(masks.shape[2]):
            mask = masks[..., region_number]
            # Write to disk
            h5file_masks.create_carray(
                "/Z_{}/{}/{}".format(z_level, sc, region_set_number),
                "region_{}".format(region_number),
                obj=mask,
                createparents=True,
                chunkshape=(128, 128),
                filters=tables.Filters(complevel=1, complib="zlib"),
            )

            # Compute the properties
            properties = skimage.measure.regionprops(mask)

            # Once a given trench has a segmentation mask, then write the properties of the masked region to the table.
            # NOTE / TODO: future versions of skimage will allow spitting out a properties object all at once, rather than lazily calculating them one at a time.
            for p in properties:
                write_properties_to_table(
                    name,
                    fov,
                    frame,
                    z_level,
                    region_set_number,
                    region_number,
                    p,
                    sc,
                    seg_params[sc],
                    row,
                    stack,
                    ch_to_img,
                )


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
        ch = ch.decode("utf-8")
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
    zeroth_channel = channels[0].decode("utf-8")
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
        ch = ch.decode("utf-8")
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
    # TODO: verify that the channels specified in the params match the available channels in the files?
    params = read_params_file(params_file)

    # HDF5 file with images & metadata
    h5file = tables.open_file(in_file, mode="r")

    # Loop the nodes immediately under Images to get the file names
    file_names = [i._v_name for i in h5file.list_nodes("/Images")]

    # Get channels from one of the files
    # Assume that all files have the same channels: otherwise, cannot process them simultaneously!
    node = h5file.get_node("/Metadata/{}".format(file_names[0]))()
    channels = get_metadata(node)["channels"]

    # Iterate the files & images
    with Pool(processes=num_cpu) as pool:
        # Make a dict of algorithms
        algo_to_func = {
            "threshold": algorithms.run_single_threshold,
            #'dual_threshold'    : algorithms.run_dual_thresholding,
            "niblack": algorithms.run_niblack_segmentation,
            "niblack_phase": algorithms.run_niblack_phase_segmentation,
            #'niblack_phase_gpu' : algorithms.run_segmentation_GPU
        }

        # Input a channel, output an algorithm
        algo_dict = {}
        for ch, p in params.items():
            algo_dict[ch] = algo_to_func[p["algorithm"]]

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

        # Analyze the images & write the masks
        print("Computing masks & measuring properties...")

        pbar = tqdm(total=total, desc="Channel stack")

        def update_pbar(*a):
            pbar.update()

        # Regions?
        if regions_file:
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
                                params,
                                channels,
                                out_dir_masks,
                                out_dir_tables,
                                algo_dict,
                                file_names,
                                regions_file,
                            ]

                            pool.apply_async(
                                run_segmentation_analysis_regions,
                                func_args,
                                callback=update_pbar,
                            )
                            # run_segmentation_analysis_regions(*func_args) # DEBUG

        else:
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
                                params,
                                channels,
                                out_dir_masks,
                                out_dir_tables,
                                algo_dict,
                                file_names,
                            ]

                            pool.apply_async(
                                run_segmentation_analysis_no_regions,
                                func_args,
                                callback=update_pbar,
                            )
                            # run_segmentation_analysis_no_regions(*func_args) # DEBUG

        pool.close()
        pool.join()

        pbar.close()

    h5file.close()

    print("Done computing masks & measuring properties.")

    # Link all the individual masks & properties files into respective H5 file, to make it easier to iterate them
    print("Linking masks...")
    link_files(out_dir_masks, "masks")
    link_files(out_dir_tables, "tables")

    # Merge the tables
    print("Merging tables...")
    in_file = os.path.join(out_dir, "TABLES/tables.h5")
    out_file = os.path.join(out_dir, "TABLES/tables_merged.h5")
    merge_tables(in_file, out_file, channels, params.keys(), file_names)

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
