#!/usr/bin/env python3

import numpy
import tables
import os
import time
from multiprocessing import Pool
from metadata import get_metadata
from tqdm import tqdm

"""
Write kymographs to disk
"""


def make_kymograph(filename, fov_number, z_level, channel, frames, h5file, trench):
    """
    Returns a 2D numpy array with all of the kymograph timepoints written sequentially, from "left" to "right".
    """
    # TODO make sure that this is the correct co-ordinate system!!
    width = trench[2] - trench[0]
    height = trench[3] - trench[1]

    # FIXME best way to figure out the dtype? assume it's float32?
    kymograph_array = numpy.empty(
        shape=(height, width * num_frames), dtype=numpy.float32
    )

    for f in frames:
        # FIXME make sure this is still correct, with F-ordered arrays
        node_path = "/Images/File_{}/FOV_{}/Frame_{}/Z_{}/{}".format(
            filename, fov_number, f, z_level, channel
        )
        node = h5file.get_node(node_path)
        kymograph_array[:, f * width : (f + 1) * width] = node[
            trench[2] : trench[3], trench[0] : trench[1]
        ]
        # TODO make sure that this is the correct co-ordinate system!!

    return kymograph_array


def query_regions_from_file(file, fov, z_level, regions_file):
    """
    Read all regions matching a given query
    Return a list of regions, and a list of trench numbers
    """
    h5file_regions = tables.open_file(regions_file, "r")
    regions_table = h5file_regions.get_node("/", "trench_coords")
    query = "(info_file == file) & (info_fov == fov) & (info_z_level == z_level)"
    search = regions_table.read_where(query)
    df = pandas.DataFrame(search)
    h5file_regions.close()

    regions = []
    trench_numbers = []
    for r in df.itertuples():
        regions.append([r.min_row, r.min_col, r.max_row, r.max_col])
        trench_numbers.append(r.info_trench_number)

    return (regions, trench_numbers)


def query_regions_from_file(
    trench_number, row_number, file, fov, z_level, regions_file
):
    """
    Read all regions matching a given query
    Return a list of regions, and a list of trench numbers
    """
    h5file_regions = tables.open_file(regions_file, "r")
    regions_table = h5file_regions.get_node("/", "trench_coords")
    query = "(info_trench_number == trench_number) \
           & (info_row_number == row_number) \
           & (info_file == file) \
           & (info_fov == fov) \
           & (info_z_level == z_level)"

    search = regions_table.read_where(query)
    df = pandas.DataFrame(search)
    h5file_regions.close()

    regions = []
    frame_numbers = []
    for r in df.itertuples():
        regions.append([r.min_row, r.min_col, r.max_row, r.max_col])
        frame_numbers.append(r.info_frame)

    return (regions, frame_numbers)


def make_kymograph_new(
    trench_number, row_number, file, fov, z_level, channel, regions_file, images_file
):
    """
    Returns a 2D numpy array with all of the kymograph timepoints written sequentially, from "left" to "right".
    """
    (trenches, frame_numbers) = query_regions_from_file(
        trench_number, row_number, file, fov, z_level, regions_file
    )

    first_trench = trenches[0]
    width = first_trench[2] - first_trench[0]
    height = first_trench[3] - first_trench[1]

    kymograph_array = numpy.empty(
        shape=(width * len(trenches), height), dtype=numpy.float32
    )

    h5file_images = tables.open_file(images_file, "r")

    # Trenches across frames (over time)
    for frame, trench in zip(frame_numbers, trenches):
        node_path = "/Images/{}/FOV_{}/Frame_{}/Z_{}/{}".format(
            file, fov, frame, z_level, channel
        )
        node = h5file_images.get_node(node_path)
        img_slice = node[trench[0] : trench[2], trench[1] : trench[3]]
        kymograph_array[frame * width : (frame + 1) * width, 0:height] = img_slice

    h5file_images.close()

    return kymograph_array


def run_kymographs(
    h5file_in, out_dir, filename, fov, channel, frames, z_level, regions
):
    """
    Input an HDF5 file and some regions information, write kymographs to a new HDF5 file.
    """

    # Open the HDF5 file
    h5file = tables.open_file(h5file_in, mode="r")

    # Create the new file for storing this FOV's kymographs
    # FIXME nested directory structure
    out_path = os.path.join(kymo_dir, "{}.h5".format(fov))
    h5file_kym = tables.open_file(out_path, mode="w")

    for i, r in enumerate(regions):
        for c in channels:
            kymograph_array = make_kymograph(
                filename, fov, z_level, c, frames, fov_h5file, r
            )

            # Foreach trench / channel, make the kymograph array & write to new file
            h5file_kym.create_carray(
                "/Trench_{}".format(i),
                c,
                obj=kymograph_array,
                chunkshape=(128, 128),
                filters=tables.Filters(complevel=1, complib="zlib"),
                createparents=True,
            )

    h5file_kym.close()
    h5file.close()


def link_files(out_dir):
    """
    Iterate the files which were just created and externally link them into a common, super h5file.
    This file stores external links to all the files stored in the HDF5 directory.
    """
    out_file = os.path.join(out_dir, "kymographs.h5")
    h5file = tables.open_file(out_file, mode="w")
    kym_dir = os.path.join(out_dir, title)

    # FIXME need to recurse through all the levels, as appropriate!
    for filename in os.listdir(kym_dir):
        (file_root, extension) = os.path.splitext(os.path.basename(filename))
        # Make an external link from this file (FOV) into the main HDF5 file
        h5file.create_external_link(
            "/",
            "{}".format(file_root),
            "kymographs/{}:/".format(filename),
            createparents=True,
        )

    h5file.close()


def main_kymographs_function(out_dir, in_file, num_cpu, regions_file):
    """ """
    # HDF5 file with images & metadata
    h5file = tables.open_file(in_file, mode="r")

    # Loop the nodes immediately under Images to get the file names
    file_names = [i._v_name for i in h5file.list_nodes("/Images")]

    total = 0
    metadata = {}
    for f in file_names:
        # Get metadata for each file
        n = h5file.get_node("/Metadata/{}".format(f))()
        m = get_metadata(n)

        # Total # of frames to be processed, for progress bar
        total += m["fields_of_view"] * m["z_levels"] * m["channels"]

        # store in a dictionary
        metadata[f] = m

    h5file.close()

    # For writing the kymographs
    os.makedirs(out_dir, exist_ok=False)

    print("Generating kymographs...")
    start = time.time()

    pbar = tqdm(total=total, desc="Kymograph")

    def update_pbar(*a):
        pbar.update()

    # Run in parallel
    with Pool(processes=num_cpu) as p:
        for filename in os.listdir(fov_dir):
            for fov in metadata[filename]["fields_of_view"]:
                regions = get_regions(regions_file, filename, fov)

                for z in metadata[filename]["z_levels"]:
                    for c in metadata[filename]["channels"]:
                        # kwargs = [h5file_in=h5file, out_dir=out_dir, filename=filename, fov=fov, channel=c, frames=metadata[filename]['frames'], z_level=z, regions=regions]
                        args = [
                            h5file,
                            out_dir,
                            filename,
                            fov,
                            c,
                            metadata[filename]["frames"],
                            z,
                            regions,
                        ]
                        p.apply_async(run_kymographs, args, callback=update_pbar)

        p.close()
        p.join()

    pbar.close()

    ### Done looping
    end = time.time()
    print(end - start)
    print("Done generating kymographs.")

    # Link in all the newly-created h5 files to the parent file
    print("Linking files...")
    link_files(out_dir)
    print("Done linking files.")
