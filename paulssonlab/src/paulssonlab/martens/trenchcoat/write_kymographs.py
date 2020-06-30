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


def make_kymograph(filename, fov_number, z_level, channel, frames, h5file, r):
    """
    Returns a 2D numpy array with all of the kymograph timepoints written sequentially, from "left" to "right".
    """
    # TODO make sure that this is the correct co-ordinate system!!
    width = r[2] - r[0]
    height = r[3] - r[1]

    # FIXME best way to figure out the dtype? assume it's float32?
    kymograph_array = numpy.empty(
        shape=(height, width * num_frames), dtype=numpy.float32
    )

    for f in frames:
        # FIXME make sure this is still correct, with F-ordered arrays
        kymograph_array[:, f * width : (f + 1) * width] = h5file.get_node(
            "/Images/File_{}/FOV_{}/Frame_{}/Z_{}/{}".format(
                filename, fov_number, f, z_level, channel
            )
        )[r[2] : r[3], r[0] : r[1]]
        # TODO make sure that this is the correct co-ordinate system!!

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
