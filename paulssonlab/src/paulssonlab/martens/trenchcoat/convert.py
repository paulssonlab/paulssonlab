#!/usr/bin/env python3

# Using nd2reader version 2.1.3
import nd2reader
import os
import tables
import time
import argparse
from multiprocessing import Pool
import pathlib

###
# NOTE: This code ignores z-level information because I never use it.
# TODO: In the future, could also make it convert a directory of ND2 files?
# TODO: pass in a range of FOV positions?
#       (maybe the user only wants a sub-set of the data, and could save time & storage)
# TODO: if min or max frames are not specified, then must provide sane defaults (i.e. 0 to max)
# TODO: it might make more sense to store the overall ND2 data in individual 1x1 arrays. More straightforward than a table.
#       This could also make it easier to dynamically add the channel names: a sub-hierarchy of the names,
#       pointing to an integer numbering (0, 1, 2...).
###

# Input an ND2 file path, the number of FOVs, a directory to write HDF5 files to, and the num. CPUs to use
# Calls the conversion function across all FOV positions in the ND2 file.
# We do not pass in an nd2reader object itself because this causes issues when dispatching to processes.
def convert_files(nd2_file_path, fov_list, hdf_dir, num_cpu, min_frame, max_frame):
    # Make a sub-directory for storing the H5 files for each FOV
    fov_path = os.path.join(hdf_dir, "FOV")
    pathlib.Path(fov_path).mkdir(parents=True, exist_ok=True)

    if num_cpu > 1:
        # https://stackoverflow.com/questions/10212445/python-map-list-item-to-function-with-arguments
        args = [
            (fov_number, fov_path, nd2_file_path, min_frame, max_frame)
            for fov_number in fov_list
        ]
        with Pool(processes=num_cpu) as p:
            p.starmap(convert_by_fov_number, args)

    else:
        for fov_number in fov_list:
            convert_by_fov_number(
                fov_number, fov_path, nd2_file_path, min_frame, max_frame
            )


# The main conversion process.
# Input the FOV number, the HDF5 directory, and the ND2 reader,
# and minimum, maximum frames to convert (and everything in between).
# Create a new HDF5 file for this FOV position (and all time frames therein).
# Copy in all the images with the compression and chunking defined below.
# FIXME is there a performance penalty in always re-initializing filters & chunk_dimensions,
# vs. init. once & passing them in? Does it matter?
def convert_by_fov_number(fov_number, fov_path, nd2_file_path, min_frame, max_frame):
    # Complevel 9 vs 1 is not worth. Very minimal gains in compression amount.
    filters = tables.Filters(complevel=1, complib="zlib")
    chunk_dimensions = (128, 128)

    reader = nd2reader.Nd2(nd2_file_path)

    # Make new HDF5 file for this fov, containing a stack of images for the different channels,
    # across specified time frames.
    out_file_new = os.path.join(fov_path, "FOV_{}.h5".format(fov_number))
    h5_new = tables.open_file(out_file_new, mode="w", title="{}".format(fov_number))

    for frame_number in range(min_frame, max_frame):
        for c in reader.channels:
            # Get the image from the ND2 file
            # Assume z_level always equals zero
            img = reader.get_image(
                channel_name=c,
                field_of_view=fov_number,
                frame_number=frame_number,
                z_level=0,
            )

            # Copy image into HDF5 file hierarchy
            h5_new.create_carray(
                "/images/Frame_{}".format(img.frame_number),
                img.channel,
                obj=img,
                title="{} raw image".format(img.channel),
                chunkshape=chunk_dimensions,
                filters=filters,
                createparents=True,
            )

    h5_new.close()
    print("FOV {} converted".format(fov_number))


# Iterate the files which were just created and externally link them into a common, super h5file.
# This file stores external links to all the files stored in the HDF5 directory.
def link_files(hdf_dir, title):
    out_file = os.path.join(hdf_dir, "{}.h5".format(title))
    h5file = tables.open_file(out_file, mode="w", title=title)

    fov_dir = os.path.join(hdf_dir, "FOV")

    for filename in os.listdir(fov_dir):
        if filename.endswith("h5"):
            (file_root, extension) = os.path.splitext(os.path.basename(filename))
            # Make an external link from this file (FOV) into the main HDF5 file
            h5file.create_external_link(
                "/",
                "{}".format(file_root),
                "FOV/{}:/".format(filename),
                createparents=True,
            )

    h5file.close()


# Quickly iterates the ND2 reader (without reading in image data) and write the metadata to 2 PyTables tables:
# 1. Global ND2 metadata
# 2. Metadata for every FOV / time frame
# TODO: apply indexing to the fov & frame number columns
# TODO: use the enumerable type for the channel names
def copy_metadata(reader, hdf_dir, title, min_frame, max_frame):
    out_file = os.path.join(hdf_dir, "{}.h5".format(title))
    h5file = tables.open_file(out_file, mode="r+", title=title)

    FOV_Frame_Time = make_fov_metadata_table_info_type()
    ND2_Metadata = make_nd2_metadata_table_info_type(reader.channels)

    # Create a table, in a separate file, for storing the timestamp information (in seconds) for each FOV / time frame
    # Then, externally link this table into the main h5file.
    # Do the same for the ND2 Reader metadata.
    metadata_table_file = os.path.join(hdf_dir, "{}_metadata.h5".format(title))
    metadata_table_h5 = tables.open_file(
        metadata_table_file, mode="w", title="Metadata"
    )

    fov_metadata_table = h5file.create_table(
        metadata_table_h5.root, "fov_metadata", FOV_Frame_Time, "FOV Metadata"
    )
    nd2_metadata_table = h5file.create_table(
        metadata_table_h5.root, "nd2_metadata", ND2_Metadata, "ND2 Metadata"
    )

    # For appending new entries
    fov_metadata_row = fov_metadata_table.row
    nd2_metadata_row = nd2_metadata_table.row

    ###

    # 1. FOV metadata
    # Zip the 4 different generators to get all 4 pieces of metadata
    num_fov = len(reader.fields_of_view)
    for i, data in enumerate(
        zip(
            reader._parser.raw_metadata.x_data,
            reader._parser.raw_metadata.y_data,
            reader._parser.raw_metadata.z_data,
            reader._parser.raw_metadata.acquisition_times,
        )
    ):

        # Convert linear index into FOV, frame
        fov = i % num_fov
        frame_number = int((i - fov) / num_fov)

        if (frame_number >= min_frame) and (frame_number <= max_frame):
            # Write the metadata into the table
            # TODO: clip or round z to 2 decimals, x & y to 1 decimal
            fov_metadata_row["info_fov"] = fov
            fov_metadata_row["info_frame"] = frame_number
            fov_metadata_row["info_x"] = data[0]
            fov_metadata_row["info_y"] = data[1]
            fov_metadata_row["info_z"] = data[2]
            fov_metadata_row["info_timestamp"] = data[3]

            fov_metadata_row.append()

    ### Done populating the metadata table
    fov_metadata_table.flush()

    # 2. ND2 metadata
    # FIXME also include the number of image cycles? (number of time frames)
    # and use max_frame - min_frame, or the full number, whichever is less
    nd2_metadata_row["info_height"] = reader.height
    nd2_metadata_row["info_width"] = reader.width
    nd2_metadata_row["info_date"] = reader.date
    nd2_metadata_row["info_fields_of_view"] = len(reader.fields_of_view)
    nd2_metadata_row["info_z_levels"] = len(reader.z_levels)
    nd2_metadata_row["info_pixel_microns"] = reader.pixel_microns

    # FIXME -- use the enumerable type
    nd2_metadata_row["info_channels"] = "\t".join(reader.channels)
    # NOTE if we set the defaults correctly, then we don't need to specify the values here
    # for c in reader.channels:
    #     nd2_metadata_row['channel_{}'.format(c)] = channels_enum[c]

    nd2_metadata_row.append()
    nd2_metadata_table.flush()

    # Link the 2 metadata tables into the super file
    h5file.create_external_link(
        "/tables",
        "fov_metadata",
        "{}_metadata.h5:/fov_metadata".format(title),
        # "{}:/fov_metadata".format(metadata_table_file),
        createparents=True,
    )

    h5file.create_external_link(
        "/tables",
        "nd2_metadata",
        "{}_metadata.h5:/nd2_metadata".format(title),
        # "{}:/nd2_metadata".format(metadata_table_file),
        createparents=True,
    )

    h5file.close()
    metadata_table_h5.close()


# Define the data types for a PyTables table.
# Table for looking up FOV, Frame,
# returning the timestamp in seconds and X, Y, Z positions.
def make_fov_metadata_table_info_type():
    # Define the PyTables column types using a dictionary
    column_types = {
        "info_fov": tables.UInt16Col(),
        "info_frame": tables.UInt16Col(),
        "info_timestamp": tables.Float32Col(),
        "info_x": tables.Float32Col(),
        "info_y": tables.Float32Col(),
        "info_z": tables.Float32Col(),
    }

    return type("FOV_Metadata", (tables.IsDescription,), column_types)


# Define the data types for a PyTables table.
# Table for storing the ND2 metadata for the whole ND2 file
# TODO: fix the channels to work with the enumerable type, rather than a fixed-length tab-delimited string.
def make_nd2_metadata_table_info_type(channel_names):
    # Define the PyTables column types using a dictionary
    column_types = {
        "info_height": tables.UInt16Col(),
        "info_width": tables.UInt16Col(),
        # TODO: make this a Time32Col, or some such? How do those work? Or is that only for time in seconds?
        # Need to look up HDF5 datetime32, datetime64 data types.
        # For now, store 19 chars: 1234567890123456789
        #                          2019-12-04 19:07:24
        "info_date": tables.StringCol(19),
        "info_fields_of_view": tables.UInt16Col(),
        "info_z_levels": tables.UInt8Col(),
        "info_pixel_microns": tables.Float32Col(),
        # FIXME for now, store channel names as tab-delimited string of length 64
        # Since the EnumCol isn't working below
        "info_channels": tables.StringCol(64),
    }

    #     # Add some column names dynamically
    #     # The idea is to add a column for each channel name
    #     # Its name is stored within the column, and its enumerated index is stored as the value.
    #     # TODO would it make more sense to have a nested column type?
    #     for c in channel_names:
    #         column_types["channel_{}".format(c)] = tables.EnumCol(channel_names, c, 'uint8'),

    #     print(column_types)
    return type("ND2_Metadata", (tables.IsDescription,), column_types)


# TODO: throw errors if the arguments don't make sense?
def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert an ND2 file to an HDF5 file equivalent."
    )
    parser.add_argument("-o", "--out-dir", type=str, help="Output directory path")
    parser.add_argument("-i", "--in-file", type=str, help="Input ND2 file path")
    parser.add_argument("-n", "--num-cpu", type=int, help="Number of CPUs to use")
    parser.add_argument("-m", "--min-frame", type=int, help="Min frame number to copy")
    parser.add_argument("-M", "--max-frame", type=int, help="Max frame number to copy")

    return parser.parse_args()


###

if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_args()

    # 1. HDF directory to store new FOV files (out-dir)
    # Store a separate HDF5 file for each FOV in this directory
    hdf_dir = args.out_dir

    # 2. Path of ND2 file to be converted (in-file)
    nd2_file_path = args.in_file

    # 3. Number of processes to run in parallel (num-cpu)
    num_cpu = args.num_cpu

    # 4. Range of time frames to be analyzed, if the user doesn't want to copy all the frames
    min_frame = args.min_frame
    max_frame = args.max_frame

    # Parse the ND2 file name to determine the title for the HDF5 file
    (nd2_file_root, extension) = os.path.splitext(os.path.basename(nd2_file_path))

    # Open the ND2 reader to get the number of FOVs, but then close it again immediately.
    reader = nd2reader.Nd2(nd2_file_path)
    fov_list = reader.fields_of_view
    ## DEBUG
    # fov_list = [0]
    reader.close()

    ### Convert the files
    print("Converting file...")
    start = time.time()
    convert_files(nd2_file_path, fov_list, hdf_dir, num_cpu, min_frame, max_frame)
    end = time.time()
    print("Done converting file, which took {} seconds.".format(end - start))

    ### Link the files
    print("Linking files...")
    link_files(hdf_dir, nd2_file_root)
    print("Done linking files...")

    ### Copy in the metadata
    # Open the ND2 reader again
    reader = nd2reader.Nd2(nd2_file_path)

    print("Copying metadata...")
    copy_metadata(reader, hdf_dir, nd2_file_root, min_frame, max_frame)
    print("Done copying metadata.")

    ### Done!
    reader.close()
