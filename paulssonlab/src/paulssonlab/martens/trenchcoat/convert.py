#!/usr/bin/env python3

# Using nd2reader version 2.1.3
import nd2reader
from multiprocessing import Pool
import pathlib
import os
import time
import tables
import numpy
import argparse
import xmltodict
import xml.etree.ElementTree as ElementTree

import warnings

# These warnings become tiresome...
warnings.simplefilter(action="ignore", category=tables.NaturalNameWarning)

###

# TODO:
# - sanity checking on input, output directories
# - alternative linking method (e.g. to encode flow lanes) (NOTE: could also make this a separate utility? Re-link afterwards.)
# - case sensitivity when checking ".nd2" -- allow ".ND2" etc.

# Main conversion process, which makes a new H5 file for each Frame
def conversion(out_dir, path, fov, frame):
    # Complevel 9 vs 1 is not worth. Very minimal gains in compression amount.
    filters = tables.Filters(complevel=1, complib="zlib")
    chunk_dimensions = (128, 128)

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
        for c in reader.channels:
            img = reader.get_image(
                channel_name=c, field_of_view=fov, frame_number=frame, z_level=z
            )

            # Copy image into HDF5 file hierarchy
            h5file.create_carray(
                "/Z_{}".format(z),
                c,
                obj=img,
                chunkshape=chunk_dimensions,
                filters=filters,
                createparents=True,
            )

    # Done!
    reader.close()
    h5file.close()


# A standard linking system, which recapitulates the directory structure
# TODO: alternative linking, for example to specify "lanes" or other types of related FOVs
def link_files(dir, frames, fields_of_view):
    out_file = os.path.join(dir, "data.h5")
    h5file_top = tables.open_file(out_file, mode="w")

    # ND2 file directories
    for nd2_entry in os.scandir(dir):
        if nd2_entry.is_dir():
            # H5 file for each ND2 file
            nd2_h5_outfile = os.path.join(dir, "{}.h5".format(nd2_entry.name))
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

            # Link metadata into the main file, too
            # Metadata will be copied in later
            h5file_top.create_external_link(
                "/Metadata",
                "File_{}".format(nd2_entry.name),
                "{}_metadata.h5:/".format(nd2_entry.name),
                createparents=True,
            )

    h5file_top.close()


# Define the data types for a PyTables table.
# Table for looking up FOV, Frame, and returning the timestamp in seconds and X, Y, Z positions.
def make_fov_metadata_table_info_type():
    # Define the PyTables column types using a dictionary
    column_types = {
        "info_fov": tables.UInt16Col(),
        "info_frame": tables.UInt16Col(),
        "info_timestamp": tables.Float64Col(),
        "info_x": tables.Float64Col(),
        "info_y": tables.Float64Col(),
        "info_z": tables.Float64Col(),
        "info_pfs_offset": tables.Float64Col(),
        "info_pfs_status": tables.BoolCol(),
    }

    return type("FOV_Metadata", (tables.IsDescription,), column_types)


# Recursively walk through a dictionary of metadata & copy over the metadata to an HDF5 file.
def dict_to_h5_metadata(dictionary, parent_node, h5file):
    for k in dictionary.keys():
        name = k.decode("utf-8")
        if name == "":
            name = "no_name"

        elem = dictionary[k]

        if type(elem) == dict:
            new_group = h5file.create_group(parent_node, name)
            dict_to_h5_metadata(elem, new_group, h5file)

        elif type(elem) == list:
            new_group = h5file.create_group(parent_node, name)
            for number, obj in enumerate(elem):
                list_to_h5_metadata(obj, new_group, h5file, "Element_{}".format(number))

        else:
            new_node = h5file.create_array(parent_node, name, obj=numpy.array(elem))


# Recursively walk through a list of metadata & copy over the metadata to an HDF5 file.
def list_to_h5_metadata(elem, parent_node, h5file, name):
    if type(elem) == dict:
        new_group = h5file.create_group(parent_node, name)
        dict_to_h5_metadata(elem, new_group, h5file)

    elif type(elem) == list:
        new_group = h5file.create_group(parent_node, name)
        for number, obj in enumerate(elem):
            list_to_h5_metadata(obj, new_group, h5file, "Element_{}".format(number))

    else:
        new_node = h5file.create_array(parent_node, name, obj=numpy.array(elem))


# After re-converting the nested OrderedDicts (JSON-like) back into XML,
# iterate all the xml elements & copy over the metadata to an HDF5 file.
# NOTE: Unicode -> ASCII conversion to handle degree symbol and micron symbol
# (problems storing unicode symbols in pytables arrays?)
def xml_to_h5_metadata_ascii(elem, parent_node, h5file, types_xml):
    for child in elem:
        # Make a new group node
        if child.attrib["runtype"] == "CLxListVariant":
            new_group = h5file.create_group(parent_node, child.tag)
            xml_to_h5_metadata_ascii(child, new_group, h5file, types_xml)

        # Make a 1x1 carray data node
        else:
            # NOTE: converted to ascii
            # To decode, use decode('unicode-escape')
            new_node = h5file.create_array(
                parent_node,
                child.tag,
                obj=numpy.array(
                    child.attrib["value"].encode("ascii", "backslashreplace"),
                    dtype=types_xml[child.attrib["runtype"]],
                ),
            )


#

# Copy ND2 metadata into an existing HDF5 hierarchy
def copy_metadata(hdf5_dir, in_file, frames, fields_of_view):
    reader = nd2reader.Nd2(in_file)
    (name, extension) = os.path.splitext(os.path.basename(in_file))

    out_file = os.path.join(hdf5_dir, "{}_metadata.h5".format(name))
    h5file = tables.open_file(out_file, mode="w")

    # Basic metadata
    copy_basic_metadata(h5file, reader, frames, fields_of_view)

    # Raw metadata
    copy_raw_metadata(h5file, reader)

    # FOV Metadata
    copy_fov_metadata(h5file, frames, fields_of_view, reader)

    h5file.close()


def copy_basic_metadata(h5file, reader, frames, fields_of_view):
    # Fields of view
    if not fields_of_view:
        fields_of_view = reader.fields_of_view

    h5file.create_array("/", "fields_of_view", obj=numpy.array(fields_of_view))

    # Z levels
    h5file.create_array("/", "z_levels", obj=numpy.array(reader.z_levels))

    # Frames
    if not frames:
        frames = reader.frames

    h5file.create_array("/", "frames", obj=numpy.array(reader.frames))

    # Height
    h5file.create_array("/", "height", obj=numpy.array(reader.height))

    # Width
    h5file.create_array("/", "width", obj=numpy.array(reader.width))

    # Pixel Microns
    h5file.create_array("/", "pixel_microns", obj=numpy.array(reader.pixel_microns))

    # Date, in UNIX time (as a float)
    h5file.create_array("/", "unix_timestamp", obj=numpy.array(reader.date.timestamp()))

    # Channels
    h5file.create_array("/", "channels", obj=numpy.array(reader.channels))


def copy_raw_metadata(h5file, reader):
    types_xml = {
        "CLxStringW": numpy.unicode,
        "lx_int32": numpy.int32,
        "lx_uint32": numpy.uint32,
        "double": numpy.float64,
        "bool": numpy.bool,
    }

    # Grabber Settings
    xml = xmltodict.unparse(reader._parser.raw_metadata.grabber_settings)
    tree = ElementTree.fromstring(xml)
    root_node = h5file.create_group(
        "/raw_metadata", "grabber_settings", createparents=True
    )
    xml_to_h5_metadata_ascii(tree, root_node, h5file, types_xml)

    # Image Attributes
    root_node = h5file.create_group(
        "/raw_metadata", "image_attributes", createparents=True
    )
    dict_to_h5_metadata(reader._parser.raw_metadata.image_attributes, root_node, h5file)

    # Image Calibration
    root_node = h5file.create_group(
        "/raw_metadata", "image_calibration", createparents=True
    )
    dict_to_h5_metadata(
        reader._parser.raw_metadata.image_calibration, root_node, h5file
    )

    # Image Metadata
    root_node = h5file.create_group(
        "/raw_metadata", "image_metadata", createparents=True
    )
    dict_to_h5_metadata(reader._parser.raw_metadata.image_metadata, root_node, h5file)

    # Image Metadata Sequence
    root_node = h5file.create_group(
        "/raw_metadata", "image_metadata_sequence", createparents=True
    )
    dict_to_h5_metadata(
        reader._parser.raw_metadata.image_metadata_sequence, root_node, h5file
    )

    # Image Text Info
    root_node = h5file.create_group(
        "/raw_metadata", "image_text_info", createparents=True
    )
    dict_to_h5_metadata(reader._parser.raw_metadata.image_text_info, root_node, h5file)

    # LUT Data
    xml = xmltodict.unparse(reader._parser.raw_metadata.lut_data)
    tree = ElementTree.fromstring(xml)
    root_node = h5file.create_group("/raw_metadata", "lut_data", createparents=True)
    xml_to_h5_metadata_ascii(tree, root_node, h5file, types_xml)

    # Custom Data
    xml = xmltodict.unparse(reader._parser.raw_metadata.custom_data)
    tree = ElementTree.fromstring(xml)
    root_node = h5file.create_group("/raw_metadata", "custom_data", createparents=True)
    xml_to_h5_metadata_ascii(tree, root_node, h5file, types_xml)


# Make a table of X, Y, Z, Timestamp, PFS information for each FOV
def copy_fov_metadata(h5file, frames, fields_of_view, reader):
    # Check whether the user specified lists of frames & fields of view to process
    # (thereby excluding those not in the list)
    if not frames:
        frames = reader.frames

    if not fields_of_view:
        fields_of_view = reader.fields_of_view

    # FOV metadata: X, Y, Z, timestamp
    FOV_Frame_Time = make_fov_metadata_table_info_type()
    fov_metadata_table = h5file.create_table(
        "/", "fov_metadata", FOV_Frame_Time, "FOV Metadata"
    )
    fov_metadata_row = fov_metadata_table.row

    num_fov = len(reader.fields_of_view)
    for i, data in enumerate(
        zip(
            reader._parser.raw_metadata.x_data,
            reader._parser.raw_metadata.y_data,
            reader._parser.raw_metadata.z_data,
            reader._parser.raw_metadata.acquisition_times,
            reader._parser.raw_metadata.pfs_offset,
            reader._parser.raw_metadata.pfs_status,
        )
    ):

        # Convert linear index into FOV, frame
        fov = i % num_fov
        frame_number = int((i - fov) / num_fov)

        # Write the metadata into the table
        if (frame_number in frames) and (fov in fields_of_view):
            # NOTE index reader on FOV, b/c there might be gaps.
            # If no gaps, then reader.fields_of_view[fov] == fov.
            fov_metadata_row["info_fov"] = reader.fields_of_view[fov]
            fov_metadata_row["info_frame"] = frame_number
            fov_metadata_row["info_x"] = data[0]
            fov_metadata_row["info_y"] = data[1]
            fov_metadata_row["info_z"] = data[2]
            fov_metadata_row["info_timestamp"] = data[3]
            fov_metadata_row["info_pfs_offset"] = data[4]
            fov_metadata_row["info_pfs_status"] = data[5]

            fov_metadata_row.append()

    # Done populating the metadata table
    fov_metadata_table.flush()


# Modified from https://rosettacode.org/wiki/Range_expansion#Python
# Input a range, such as: 1,4-7
# Return a list with all elements within the range.
def range_expand(range_string):
    result = []

    for r in range_string.split(","):
        # Start at 1, because the first number might be negative,
        # meaning that the first char is a minus (dash) sign.
        if "-" in r[1:]:
            first, second = r[1:].split("-", 1)
            result += range(int(r[0] + first), int(second) + 1)

        else:
            result.append(int(r))

    return result


# TODO: throw errors if the arguments don't make sense?
def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert an ND2 file to an HDF5 file equivalent."
    )

    parser.add_argument("-o", "--out-dir", type=str, help="Output directory.")
    parser.add_argument(
        "-i", "--in-dir", type=str, help="Input directory, with 1 or more ND2 files."
    )
    parser.add_argument(
        "-n",
        "--num-cpu",
        type=int,
        help="Number of CPUs to use. Default: use all CPUs.",
    )
    parser.add_argument(
        "-f",
        "--frames",
        type=str,
        help="List of frames to copy. Default: copy all frames.",
    )
    parser.add_argument(
        "-F", "--fovs", type=str, help="List of FOVs to copy. Default: copy all FOVs."
    )

    return parser.parse_args()


###

if __name__ == "__main__":
    # Parse command-line arguments
    # TODO do some basic sanity checking on the inputs (not empty, etc.)
    args = parse_args()

    # 1. HDF directory to store new FOV files (out-dir)
    # Store a separate HDF5 file for each FOV in this directory
    hdf_dir = args.out_dir

    # 2. Path of ND2 file to be converted (in-file)
    nd2_dir = args.in_dir

    # 3. Number of processes to run in parallel
    if args.num_cpu:
        num_cpu = args.num_cpu
    else:
        num_cpu = None  # defaults to use all available CPUs

    # 4. Parse the frame and fov ranges
    # a. Range of time frames to be analyzed, if the user doesn't want to copy all the frames
    if args.frames:
        frames = range_expand(args.frames)
    else:
        frames = None

    # b. Range of FOVs to be analyzed
    if args.fovs:
        fields_of_view = range_expand(args.fovs)
    else:
        fields_of_view = None

    ###

    # Write the images
    print("Writing files...")
    start = time.time()

    with Pool(processes=num_cpu) as p:
        # Iterate the ND2 files in the directory
        # Use sorted order, to be more consistent.
        files = sorted([file for file in os.scandir(nd2_dir)], key=lambda f: f.path)

        for in_file in files:
            if in_file.path.endswith(".nd2"):
                nd2_file = nd2reader.Nd2(in_file.path)

                if not fields_of_view:
                    fields_of_view = nd2_file.fields_of_view

                if not frames:
                    frames = nd2_file.frames

                z_levels = nd2_file.z_levels

                nd2_file.close()

                for fov in fields_of_view:
                    for frame in frames:
                        args = [hdf_dir, in_file.path, fov, frame]
                        p.apply_async(conversion, args)

        p.close()
        p.join()

    end = time.time()
    print("Done writing images, which took {} seconds.".format(end - start))

    # Link the files
    print("Linking files...")
    start = time.time()
    link_files(hdf_dir, frames, fields_of_view)
    end = time.time()
    print("Done linking files, which took {} seconds.".format(end - start))

    # Copy the metadata
    print("Copying metadata...")
    start = time.time()

    # Don't use more processes than files
    files = []
    for f in os.scandir(nd2_dir):
        if f.path.endswith(".nd2"):
            files.append(f)

    if num_cpu:
        if num_cpu > len(files):
            num_cpu = len(files)

    with Pool(processes=num_cpu) as p:
        for in_file in files:
            args = [hdf_dir, in_file.path, frames, fields_of_view]
            p.apply_async(copy_metadata, args)

        p.close()
        p.join()

    end = time.time()
    print("Done copying metadata, which took {} seconds.".format(end - start))

    # Done!
