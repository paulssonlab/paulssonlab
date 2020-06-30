#!/usr/bin/env python3

# Using nd2reader version 2.1.3, with modifications for loading F-ordered channel stacks
import nd2reader
from multiprocessing import Pool
import pathlib
import os
import tables
import numpy
from tqdm import tqdm

import metadata

import warnings

# These warnings become tiresome...
warnings.simplefilter(action="ignore", category=tables.NaturalNameWarning)

###

# TODO:
# - sanity checking on input, output directories
# - Before copying data, check that the channels are all the same in all the files:
#   otherwise, it's inappropriate to convert these files together!
# - Could also enforce same dimensions, if this helps with browsing e.g. using Napari?
#   But it mightn't make sense to enforce FOV count, because that could conceivably vary from one file to another.


def make_nd2_list(nd2_dir):
    """
    Input a directory, output a list of nd2 file paths.
    """
    files = []
    for f in os.scandir(nd2_dir):
        if f.path.lower().endswith(".nd2"):
            files.append(f.path)

    return files


def conversion(out_dir, path, fov, frame):
    """
    Main conversion process, which makes a new H5 file for each Frame
    """
    # Create all parent directories
    (file_root, extension) = os.path.splitext(os.path.basename(path))
    dir_path = "{}/{}/FOV_{}/".format(out_dir, file_root, fov, frame)
    pathlib.Path(dir_path).mkdir(parents=True, exist_ok=True)

    # Create the Frame HDF5 file, containing all Zs and channels at that time frame
    out_file = os.path.join(dir_path, "Frame_{}.h5".format(frame))
    h5file = tables.open_file(out_file, mode="w")

    # Open the ND2 File & loop the Z/channel stack
    reader = nd2reader.Nd2(path)

    for z in reader.z_levels:
        # Read in a 3D array, width x height x channels, in F-order
        stack = reader.get_image_stack(field_of_view=fov, frame_number=frame, z_level=z)

        # NOTE this iteration will work if the channel ordering is preserved, which it should be!
        for i, c in enumerate(reader.channels):
            # Copy image into HDF5 file hierarchy
            h5file.create_carray(
                "/Z_{}".format(z),
                c,
                obj=stack[..., i],
                chunkshape=(128, 128),
                # Complevel 9 vs 1 is not worth it. Very minimal gains in compression amount.
                filters=tables.Filters(complevel=1, complib="zlib"),
                createparents=True,
            )

    # Done!
    reader.close()
    h5file.close()


def link_files(in_dir, frames, fields_of_view):
    """
    A standard linking system, which recapitulates the directory structure
    """
    out_file = os.path.join(in_dir, "data.h5")
    h5file_top = tables.open_file(out_file, mode="w")

    out_file_metadata = os.path.join(in_dir, "metadata.h5")
    h5file_metadata = tables.open_file(out_file_metadata, mode="w")

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
                    # Parse the number of the FOV & check that it's in the list of desired FOVs
                    fov_num = int(fov_entry.name.split("_")[1])
                    if fov_num in fields_of_view:
                        h5file_fov = tables.open_file(
                            "{}.h5".format(fov_entry.path),
                            mode="w",
                            title=fov_entry.name,
                        )

                        # The individual H5 files, for each time frame (containing Z/channel stacks)
                        for frame_entry in os.scandir(fov_entry):
                            if frame_entry.is_file():
                                (name, extension) = os.path.splitext(
                                    os.path.basename(frame_entry.name)
                                )
                                frame_num = int(name.split("_")[1])
                                if (frame_num in frames) and (extension == ".h5"):
                                    h5file_fov.create_external_link(
                                        "/",
                                        "{}".format(name),
                                        "{}/{}:/".format(
                                            fov_entry.name, frame_entry.name
                                        ),
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
                "File_{}".format(nd2_entry.name),
                "{}.h5:/".format(nd2_entry.name),
                createparents=True,
            )

        # Link metadata into the overall metadata file
        elif nd2_entry.name.endswith("_metadata.h5"):
            name = nd2_entry.name.split("_metadata.h5")[0]
            h5file_metadata.create_external_link(
                "/", "File_{}".format(name), "{}_metadata.h5:/".format(name)
            )

    h5file_metadata.close()

    # And the final link of the metadata into the main file
    h5file_top.create_external_link("/", "Metadata", "metadata.h5:/")

    h5file_top.close()


def main_conversion_function(hdf_dir, nd2_dir, num_cpu, frames, fields_of_view):
    """
    Copy the metadata, write the images, and link the H5 files
    """
    pathlib.Path(hdf_dir).mkdir(parents=True, exist_ok=True)

    files = make_nd2_list(nd2_dir)

    # Metadata
    pbar = tqdm(total=len(files), desc="File metadata")

    def update_pbar(*a):
        pbar.update()

    with Pool(processes=num_cpu) as p:
        for in_file in files:
            args = [hdf_dir, in_file, frames, fields_of_view]
            p.apply_async(metadata.copy_metadata, args, callback=update_pbar)

        p.close()
        p.join()

    pbar.close()

    # Write the images
    total_frames = 0

    # Iterate the ND2 files in the directory
    # Use sorted order, to be more consistent.
    files = sorted(
        [file for file in os.scandir(nd2_dir) if file.path.lower().endswith(".nd2")],
        key=lambda f: f.path,
    )

    # First, loop the files & calculate how many frames there are total (for the progress bar):
    for in_file in files:
        nd2_file = nd2reader.Nd2(in_file)

        if not fields_of_view:
            fields_of_view = nd2_file.fields_of_view

        if not frames:
            frames = nd2_file.frames

        total_frames += len(fields_of_view) * len(frames)

    pbar_frames = tqdm(total=total_frames, desc="Frame image data")

    def update_pbar_frames(*a):
        pbar_frames.update()

    # Now, loop again
    with Pool(processes=num_cpu) as p:
        for in_file in files:
            nd2_file = nd2reader.Nd2(in_file)

            if not fields_of_view:
                fields_of_view = nd2_file.fields_of_view

            if not frames:
                frames = nd2_file.frames

            nd2_file.close()

            for fov in fields_of_view:
                for frame in frames:
                    args = [hdf_dir, in_file, fov, frame]
                    p.apply_async(conversion, args, callback=update_pbar_frames)

        p.close()
        p.join()

    pbar_frames.close()

    ## Link the files
    print("Linking files...")
    link_files(hdf_dir, frames, fields_of_view)

    ## FIXME TODO implement errors here!
    # Or maybe just print a warning...
    ## Run a check that the channels in each of the files' metadata are the same,
    ## and abort conversion process if not.
    # channels_equal = metadata_channels_equal(h5file_metadata.list_nodes("/"))

    # if channels_equal:
    # pass # no error!
    # else:
    # pass # error!

    # Done!
