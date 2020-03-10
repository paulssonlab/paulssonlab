#!/usr/bin/env python3

import numpy
import tables
import os
import argparse
import time
from multiprocessing import Pool

### Write kymographs to disk
#
#


def make_kymograph(
    fov_number, channel, y_dimension, trench_width, num_frames, h5file, ranges
):
    kymograph_array = numpy.empty(
        shape=(y_dimension, trench_width * num_frames), dtype=numpy.uint16
    )

    for f in range(num_frames):
        kymograph_array[:, f * trench_width : (f + 1) * trench_width] = h5file.get_node(
            "/images/Frame_{}/{}".format(f, channel)
        )[crop_top:crop_bottom, ranges[0] : ranges[1]]

    return kymograph_array


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


def run_kymographs(
    fov_file_path, kymo_dir, crop_top, crop_bottom, y_dimension, channels, identifier
):

    # Extract the FOV number from the file name
    base_name = os.path.basename(fov_file_path)
    (root, extension) = os.path.splitext(base_name)
    fov = int(root.split("_")[1])

    print("FOV_{}".format(fov))

    # For writing the segmentation masks
    filters = tables.Filters(complevel=1, complib="zlib")
    chunk_dimensions = (128, 128)

    # Open this FOV file
    fov_h5file = tables.open_file(fov_file_path, mode="r+")

    # Open the detected trenches table
    ranges = get_trench_coordinates(
        "/tables/trench_coordinates", identifier, fov_h5file
    )
    # TODO infer automatically?
    trench_width = 30

    # How many frames?
    num_frames = 0
    for time_node in fov_h5file.root.images._f_iter_nodes():
        num_frames += 1

    # Create the new file for storing this FOV's kymographs
    out_path = os.path.join(kymo_dir, "kymographs_FOV_{}.h5".format(fov))
    h5file_kym = tables.open_file(out_path, mode="w", title="Kymographs")

    # Iterate each trench
    for i, tr in enumerate(ranges):

        # Iterate each channel
        for c in channels:
            kymograph_array = make_kymograph(
                fov, c, y_dimension, trench_width, num_frames, fov_h5file, tr
            )

            # Foreach trench / channel, make the kymograph array & write to new file
            h5file_kym.create_carray(
                "/Trench_{}".format(i),
                c,
                obj=kymograph_array,
                title="{}/{} kymograph".format(i, c),
                chunkshape=chunk_dimensions,
                filters=filters,
                createparents=True,
            )

    h5file_kym.close()


# Iterate the files which were just created and externally link them into a common, super h5file.
# This file stores external links to all the files stored in the HDF5 directory.
def link_files(hdf_dir, title):
    out_file = os.path.join(hdf_dir, "{}.h5".format(title))
    h5file = tables.open_file(out_file, mode="w", title=title)

    kym_dir = os.path.join(hdf_dir, title)

    for filename in os.listdir(kym_dir):
        if filename.endswith("h5"):
            (file_root, extension) = os.path.splitext(os.path.basename(filename))
            # Make an external link from this file (FOV) into the main HDF5 file
            h5file.create_external_link(
                "/",
                "{}".format(file_root),
                "kymographs/{}:/".format(filename),
                createparents=True,
            )

    h5file.close()


# TODO: throw errors if the arguments don't make sense?
def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate kymographs & write them to disk."
    )
    parser.add_argument("-o", "--out-file", type=str, help="Output HDF5 file")
    parser.add_argument(
        "-i",
        "--in-dir",
        type=str,
        help="Input directory containing HDF5 files & directories",
    )
    parser.add_argument("-n", "--num-cpu", type=int, help="Number of CPUs to use")
    parser.add_argument("-I", "--identifier", type=str, help="Trench row identifier")
    parser.add_argument("-T", "--crop-top", type=int, help="Crop top")
    parser.add_argument("-B", "--crop-bottom", type=int, help="Crop bottom")

    return parser.parse_args()


###

if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_args()

    in_dir = args.in_dir
    out_file = args.out_file
    num_cpu = args.num_cpu
    identifier = args.identifier
    crop_top = args.crop_top
    crop_bottom = args.crop_bottom

    y_dimension = crop_bottom - crop_top

    # TODO: read channels from metadata?
    channels = ["BF", "mCherry", "YFP", "CFP"]

    fov_dir = os.path.join(in_dir, "FOV")
    files = os.listdir(fov_dir)
    out_dir = os.path.join(in_dir, "kymographs")
    os.makedirs(out_dir, exist_ok=True)

    print("Generating kymographs...")
    start = time.time()

    # Run in parallel
    if num_cpu > 1:
        args = [
            (
                os.path.join(fov_dir, filename),
                out_dir,
                crop_top,
                crop_bottom,
                y_dimension,
                channels,
                identifier,
            )
            for filename in os.listdir(fov_dir)
        ]

        # TODO: weed out any files which do not end in h5

        with Pool(processes=num_cpu) as p:
            p.starmap(run_kymographs, args)
    else:
        for filename in files:
            if filename.endswith("h5"):
                run_kymographs(
                    os.path.join(fov_dir, filename),
                    out_dir,
                    crop_top,
                    crop_bottom,
                    y_dimension,
                    channels,
                    identifier,
                )

    ### Done looping all FOV
    end = time.time()
    print(end - start)
    print("Done generating kymographs.")

    # Link in all the newly-created h5 files to the parent file
    print("Linking files...")
    link_files(in_dir, "kymographs")
    print("Done linking files.")
