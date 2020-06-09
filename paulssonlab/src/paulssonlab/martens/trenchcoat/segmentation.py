#!/usr/bin/env python3

import os
import time
import pathlib
import tables
from multiprocessing import Pool
from tqdm import tqdm

import algorithms
import arrayfire_algorithms
from properties import write_properties_to_table, merge_tables
from params import read_params_file

### TODO:
#  1. GPU seg algorithm as one of the options. (ArrayFire + skimage watershed).
#  2. Pass in min. region size. Use pixel_microns to help convert from true area in microns^2 to region size in pixels.
#     Makes this parameter magnification-independent.

### Perform cell segmentation & measure fluorescence intensities in microscope images, with support for sub-regions (e.g. "trenches").
###

# Input a node in the HDF5 metadata section, return a dict. with the following metadata
def get_metadata(n):
    metadata = {
        "channels": n.channels.read(),
        "fields_of_view": n.fields_of_view.read(),
        "frames": n.frames.read(),
        "width": n.width.read(),
        "height": n.height.read(),
        "pixel_microns": n.pixel_microns.read(),
        "unix_timestamp": n.unix_timestamp.read(),
        "z_levels": n.z_levels.read(),
    }

    return metadata


# Link the masks or tables files into a single HDF5 file
def link_files(in_dir, file_name):
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
    # Dir containing new HDF5 files, write results to
    out_dir_masks = os.path.join(out_dir, "MASKS")
    out_dir_tables = os.path.join(out_dir, "TABLES")

    # Create the masks & tables directories
    pathlib.Path(out_dir_masks).mkdir(parents=True, exist_ok=True)
    pathlib.Path(out_dir_tables).mkdir(parents=True, exist_ok=True)

    # Segmentation parameters, in YAML
    # TODO: verify that the channels specified in the params match the available channels in the files?
    params = read_params_file(params_file)

    ###

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

        # Regions?
        if regions_file:
            func = algorithms.run_segmentation_analysis_regions
            h5file_regions_file = regions_file

        else:
            func = algorithms.run_segmentation_analysis_noregions
            h5file_regions_file = None

        # Analyze the images & write the masks
        print("Computing masks & measuring properties...")

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

        pbar = tqdm(total=total, desc="Channel stack")

        def update_pbar(*a):
            pbar.update()

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
                            h5file_regions_file,
                            file_names,
                        ]

                        pool.apply_async(func, func_args, callback=update_pbar)
                        # func(*func_args) # DEBUG

        pool.close()
        pool.join()

        pbar.close()

    h5file.close()

    ###

    # Link all the individual masks & properties files into respective H5 file, to make it easier to iterate them
    print("Linking masks...")
    start = time.time()

    link_files(out_dir_masks, "masks")
    link_files(out_dir_tables, "tables")

    end = time.time()
    print("Done linking masks ({} seconds).".format(end - start))

    ###

    ### Merge the tables
    print("Merging tables...")
    start = time.time()

    in_file = os.path.join(out_dir, "TABLES/tables.h5")
    out_file = os.path.join(out_dir, "TABLES/tables_merged.h5")
    merge_tables(in_file, out_file, channels, params.keys(), file_names)

    end = time.time()
    print("Done merging tables ({} seconds).".format(end - start))
