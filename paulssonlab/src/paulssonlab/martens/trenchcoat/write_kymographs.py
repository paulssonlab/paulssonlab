#!/usr/bin/env python

import os
import tables
import numpy
import pandas
import pathlib
from multiprocessing import Pool
from tqdm import tqdm
from params import read_params_string
from napari_browse_hdf5 import (
    get_largest_extents_hdf5,
    metadata_attributes_equal,
    metadata_array_equal,
)
from lineages import relabel_mask_complete

"""
Write kymographs to disk. This only includes the original intensity images, not mask, regions or other "layers."

The goal of this functionality is to speed up visualization in e.g. Napari. I already have code that generates kymographs on-the-fly,
but it can be very slow (possibly limited by image transfer speeds over the network).

The intended use is *not* to use the kymographs for measurements. The philosophy of TrenchCoat is to always
refer back to the original intensity images for extracting information.

Thus pre-rendering these kymographs can speed up visualization time, at the cost of generating the kymographs up front
& having to store them on disk somewhere.

Requires first correcting for stage drift (see: renumber_trenches.py, trenchcoat correct_trench_numbering).

TODO: other things for masks

- MASKS/Parameters contains seg_params.yaml, which includes the names of the seg channels
- when calling, specify whether is a masks or an intensity images kymograph. changes:
    - add the option to the main args
    - check this arg
    - whether to use channels, or seg_channels; & which to write to the metadata
    - call write_this_group, or write_masks
- Should the chunk size be set to exactly equal a trench? This way, grabbing a trench will always be exactly one chunk!
"""


def main_kymographs_function(images_file, masks_file, regions_file, out_dir, num_cpu):
    """
    Main function, invoked from the command-line.
    """
    # Load the _entire_ set of regions (trenches) into a dataframe, and pass it to each process.
    # Call it df_regions. During the loop, a subset of the dataframe will be pulled out, gradually
    # narrowing down to the final set.
    h5file_reg = tables.open_file(regions_file, "r")
    table = h5file_reg.get_node("/regions_table")
    df_regions = pandas.DataFrame(table.read())
    h5file_reg.close()

    # Convert bytes to str for the "info_file" column for easier string comparisons
    df_regions["info_file"] = df_regions["info_file"].map(lambda x: x.decode("utf-8"))

    # File with images & metadata
    h5file_img = tables.open_file(images_file, "r")

    # Get the largest extents across all the nd2 files within the h5 file
    attributes = ["fields_of_view", "frames", "z_levels"]
    extents = {}
    for a in attributes:
        (smallest, largest) = get_largest_extents_hdf5(h5file_img, a)
        if smallest == largest:
            extents[a] = [smallest]
        else:
            extents[a] = [i for i in range(smallest, largest + 1)]

    # Define the channels
    channels = metadata_array_equal(h5file_img, "channels")
    channels = [c.decode("utf-8") for c in channels]

    # Get the names of the "File" nodes
    file_nodes = [x._v_name for x in h5file_img.list_nodes("/Images")]

    # Done reading metadata from the images HDF5 file
    # This file will be re-opened later, to read the images
    h5file_img.close()

    # Process kymographs for segmentation masks? If so, then grab the seg channel names
    if masks_file:
        # Open the file
        h5file_masks = tables.open_file(masks_file, "r")

        # Load the yaml file with the segmentation parameters
        seg_params_node = h5file_masks.get_node("/Parameters", "seg_params.yaml")
        seg_params_str = seg_params_node.read().tobytes().decode("utf-8")
        seg_params_dict = read_params_string(seg_params_str)
        seg_channels = [k for k in seg_params_dict.keys()]

        # Done for now
        h5file_masks.close()

    # Total number of frames
    num_frames = len(extents["frames"])

    # Get a sample trench to determine the width & height
    # Assume all trenches have equal dimensions
    sample_trench = df_regions.iloc[0]
    tr_width = sample_trench.max_row - sample_trench.min_row
    height = sample_trench.max_col - sample_trench.min_col

    # The width of each kymograph image
    kymo_width = tr_width * num_frames

    # Set up a progress bar
    # NOTE Assumes same num. fovs, etc. in every file node!
    total = len(file_nodes) * len(extents["fields_of_view"]) * len(extents["z_levels"])
    pbar = tqdm(total=total, desc="FOV #")

    def update_pbar(*a):
        pbar.update()

    # with Pool(processes=num_cpu) as pool:
    # for file in file_nodes:
    ## Narrow down the regions dataframe to just this file
    # df_reg_f = df_regions[ df_regions["info_file"] == file]

    # for fov in extents["fields_of_view"]:
    ## Narrow down the regions dataframe to just this fov
    # df_reg_fov = df_reg_f[ df_reg_f["info_fov"] == fov]

    # for z_level in extents["z_levels"]:
    ## Narrow down the regions dataframe to just this z_level
    # df_reg_z = df_reg_fov[ df_reg_fov["info_z_level"] == z_level]

    ## Masks
    # if masks_file:
    # args = [masks_file, out_dir, file, fov, z_level, seg_channels, extents["frames"], df_reg_z, kymo_width, tr_width, height]
    # pool.apply_async(write_masks, args, callback=update_pbar)
    ##write_masks(*args) # DEBUG

    ## Intensity images
    # else:
    # args = [images_file, out_dir, file, fov, z_level, channels, extents["frames"], df_reg_z, kymo_width, tr_width, height]
    # pool.apply_async(write_this_group, args, callback=update_pbar)
    ##write_this_group(*args) # DEBUG

    # pool.close()
    # pool.join()

    # Now, link up all the sub-HDF5 files into a parent file.
    link_files(out_dir, masks_file)

    # Write the seg params, or some metadata
    if masks_file:
        write_seg_params(out_dir, seg_params_str)
    else:
        write_metadata(
            out_dir,
            extents["fields_of_view"],
            extents["z_levels"],
            extents["frames"],
            channels,
            tr_width,
            height,
        )

    # Done!


def write_seg_params(out_dir, seg_params_str):
    """
    Write a copy of the segmentation parameters to disk.
    """
    h5file_masks = tables.open_file(os.path.join(out_dir, "kymographs.h5"), "r+")
    h5file_masks.create_array(
        "/Parameters",
        "seg_params.yaml",
        obj=numpy.array(seg_params_str),
        createparents=True,
    )
    h5file_masks.close()


def link_files(out_dir, masks_file):
    """
    After writing all the kymograph intensity images to separate sub-HDF5 files,
    use HDF5 internal links to link them up into a single hierarchy.
    """
    print("Linking files...")
    out_file = os.path.join(out_dir, "kymographs.h5")
    h5file_top = tables.open_file(out_file, mode="w")

    # ND2 file directories
    for nd2_entry in os.scandir(out_dir):
        if nd2_entry.is_dir():
            # H5 file for each ND2 file
            nd2_h5_outfile = os.path.join(out_dir, "{}.h5".format(nd2_entry.name))
            h5file_nd2file = tables.open_file(
                nd2_h5_outfile, mode="w", title=nd2_entry.name
            )

            # FOV directories
            for fov_entry in os.scandir(nd2_entry):
                if fov_entry.is_dir():
                    # Parse the number of the FOV
                    fov_num = int(fov_entry.name.split("_")[1])

                    h5file_fov = tables.open_file(
                        "{}.h5".format(fov_entry.path), mode="w", title=fov_entry.name
                    )

                    # The individual H5 files, for each time frame (containing Z/Row/Region/channel)
                    for frame_entry in os.scandir(fov_entry):
                        if frame_entry.is_file():
                            (name, extension) = os.path.splitext(
                                os.path.basename(frame_entry.name)
                            )
                            frame_num = int(name.split("_")[1])
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
                "/Images",
                "{}".format(nd2_entry.name),
                "{}.h5:/".format(nd2_entry.name),
                createparents=True,
            )

    if not masks_file:
        # And the final link of the metadata into the main file
        h5file_top.create_external_link("/", "Metadata", "metadata.h5:/")

    # Done!
    h5file_top.close()
    print("Done linking files.")


def write_metadata(
    out_dir, fields_of_view, z_levels, frames, channels, tr_width, height
):
    """
    Write some useful metadata for the on-disk kymographs:
    - height
    - width of entire kymograph
    - channels
    - fields of view
    - z levels
    - trenches for each fov/z_level
    (store as a what? bunch of lists, each one in its own node? or just a number, assuming it starts at zero?
    (which mightn't be true if they are corrected for drift))

    TODO Make it nested, for each file, just like regular datasets?
    Good way to copy all the original metadata?
    3 data which are different:
    - frames
    - height
    - width

    And 2 data which are new:
    - trench rows (for each FOV)
    - trench numbers (for each FOV)
    """
    print("Writing metadata...")
    out_path = os.path.join(out_dir, "metadata.h5")
    h5file = tables.open_file(out_path, "w")

    # Fields of view
    h5file.create_array("/", "fields_of_view", obj=numpy.array(fields_of_view))

    # Z levels
    h5file.create_array("/", "z_levels", obj=numpy.array(z_levels))

    # Frames -- FIXME does this have any meaning for a kymograph???
    h5file.create_array("/", "frames", obj=numpy.array(frames))

    # Height -- NOTE different from source file (is trench height)
    h5file.create_array("/", "height", obj=numpy.array(height))

    # Width -- NOTE different from source file (is trench width)
    h5file.create_array("/", "width", obj=numpy.array(tr_width))

    # Channels
    h5file.create_array("/", "channels", obj=numpy.array(channels))

    # TODO write out arrays of trench rows & numbers for each fov/z??

    # NOTE decided against this strategy
    ## is_kymograph ~ special flag to distinguish from ordinary HDF5 image datasets.
    # h5file.create_array("/", "is_kymograph", obj=numpy.array(True))

    h5file.close()
    print("Done writing metadata!")


def write_this_group(
    images_file,
    out_dir,
    file,
    fov,
    z_level,
    channels,
    frames,
    df_reg_z,
    kymo_width,
    tr_width,
    height,
):
    """
    The function which takes care of writing the various channels as separate nodes to a new HDF5 file.
    Only processes a single region (trench).
    Use df_reg_z dataframe to make it easier to iterate over the rows & regions.
    Pass in number of frames, kymo_width, tr_width, height to set up the final image size.
    """
    # Write the results to this file
    dir_path = "{}/{}/FOV_{}/".format(out_dir, file, fov)
    pathlib.Path(dir_path).mkdir(parents=True, exist_ok=True)
    out_path = os.path.join(dir_path, "Z_{}.h5".format(z_level))
    h5file_out = tables.open_file(out_path, "w")

    # Open the HDF5 file with the images
    h5file_img = tables.open_file(images_file, "r")

    # Look up in the table to determine how many trench rows in this FOV
    trench_rows = df_reg_z["info_row_number"].unique()
    max_trenches = df_reg_z["corrected_trench_label"].max() + 1

    for ch in channels:
        # Initialize a single ndarray. index by row-trench-[x,y] <in kymograph coordinates>.
        # Load intensity image, load trench coords, & copy into the ndarray.
        # Then, loop over trench rows & trenches, & write out the contents.
        # For now, assume that trench_rows always starts at 0 and increases by 1.
        # So the indexing will fit neatly.
        framebuffer = numpy.zeros(
            shape=(len(trench_rows), max_trenches, kymo_width, height),
            dtype=numpy.uint16,
        )

        # Enumerate: index gives position, while fr gives identity (true numbering).
        # Do this in case some frames were cropped out during ND2 -> HDF5 conversion.
        for index, fr in enumerate(frames):
            # Load the intensity image
            node_path = "/Images/{}/FOV_{}/Frame_{}/Z_{}/{}".format(
                file, fov, fr, z_level, ch
            )
            image = h5file_img.get_node(node_path)

            df_reg_frame = df_reg_z[df_reg_z["info_frame"] == fr]

            for row in trench_rows:
                # Load all the regions for this particular frame
                # Parse df_reg_z to get the trench bbox coords for this file-fov-frame-z-row
                # The input dataframe was filtered to this file-fov-z.
                df_reg_row = df_reg_frame[df_reg_frame["info_row_number"] == row]

                for r in df_reg_row.itertuples():
                    # Read the slice out of the intensity image, using the bbox coords on the trench
                    # Write the contents into the framebuffer
                    framebuffer[
                        row,
                        r.corrected_trench_label,
                        index * tr_width : (index + 1) * tr_width,
                        0:height,
                    ] = image[
                        r.min_row + r.info_x : r.max_row + r.info_x,
                        r.min_col - r.info_y : r.max_col - r.info_y,
                    ]
                    # NOTE take into account the drift correction back to the raw image

        # Write the new images to disk
        for row in range(framebuffer.shape[0]):
            for tr in range(framebuffer.shape[1]):
                h5file_out.create_carray(
                    "/row_{}/region_{}".format(row, tr),
                    ch,
                    obj=framebuffer[row, tr, :],
                    chunkshape=(128, 128),
                    # Complevel 9 vs 1 is not worth it. Very minimal gains in compression amount.
                    filters=tables.Filters(complevel=1, complib="zlib"),
                    createparents=True,
                )

    # Close up shop
    h5file_img.close()
    h5file_out.close()


def write_masks(
    masks_file,
    out_dir,
    file,
    fov,
    z_level,
    seg_channels,
    frames,
    df_reg_z,
    kymo_width,
    tr_width,
    height,
):
    """
    Take labeled trenches (from a cell segmentation pipeline) & write their kymographs.
    Unlike intensity images, each labeled trench is stored in its own, small node in an HDF5 file.
    Use df_reg_z dataframe to make it easier to iterate over the rows & regions (rather than iterating
    through HDF5 file nodes, which is surprisingly hairy).

    Pass in number of frames, kymo_width, tr_width, height to set up the final image size.

    TODO there is a lot of overlap with the intensity images version. Is there a good way
    to merge the two functions?
    """
    # Write the results to this file
    dir_path = "{}/{}/FOV_{}/".format(out_dir, file, fov)
    pathlib.Path(dir_path).mkdir(parents=True, exist_ok=True)
    out_path = os.path.join(dir_path, "Z_{}.h5".format(z_level))
    h5file_out = tables.open_file(out_path, "w")

    # Open the HDF5 file with the masks
    h5file_masks = tables.open_file(masks_file, "r")

    # Look up in the table to determine how many trench rows in this FOV
    trench_rows = df_reg_z["info_row_number"].unique()
    max_trenches = df_reg_z["corrected_trench_label"].max() + 1

    for sc in seg_channels:
        # Initialize a single ndarray. index by row-trench-[x,y] <in kymograph coordinates>.
        # Load mask, load trench coords, & copy into the ndarray.
        # Then, loop over trench rows & trenches, & write out the contents.
        # For now, assume that trench_rows always starts at 0 and increases by 1.
        # So the indexing will fit neatly.
        framebuffer = numpy.zeros(
            shape=(len(trench_rows), max_trenches, kymo_width, height),
            dtype=numpy.uint16,
        )

        # Enumerate: index gives position, while fr gives identity (true numbering).
        # Do this in case some frames were cropped out during ND2 -> HDF5 conversion.
        for index, fr in enumerate(frames):
            df_reg_frame = df_reg_z[df_reg_z["info_frame"] == fr]

            for row in trench_rows:
                # Load all the regions for this particular frame
                # Parse df_reg_z to get the trench bbox coords for this file-fov-frame-z-row
                # The input dataframe was filtered to this file-fov-z.
                df_reg_row = df_reg_frame[df_reg_frame["info_row_number"] == row]

                for r in df_reg_row.itertuples():
                    # Read the mask
                    node_path = "/masks/{}/FOV_{}/Frame_{}/Z_{}/{}/{}/region_{}".format(
                        file, fov, fr, z_level, sc, row, r.info_trench_number
                    )
                    mask = h5file_masks.get_node(node_path).read()

                    # Write the contents into the framebuffer
                    framebuffer[
                        row,
                        r.corrected_trench_label,
                        index * tr_width : (index + 1) * tr_width,
                        0:height,
                    ] = mask

        # Write the new images to disk
        for row in range(framebuffer.shape[0]):
            for tr in range(framebuffer.shape[1]):
                h5file_out.create_carray(
                    "/row_{}/region_{}".format(row, tr),
                    sc,
                    obj=framebuffer[row, tr, :],
                    # TODO make the chunk shape match the trench dimensions?
                    # Which is width, & height?
                    chunkshape=(128, 128),
                    # Complevel 9 vs 1 is not worth it. Very minimal gains in compression amount.
                    filters=tables.Filters(complevel=1, complib="zlib"),
                    createparents=True,
                )

    # Close up shop
    h5file_masks.close()
    h5file_out.close()
