#!/usr/bin/env python3
import os

import numpy
import tables
from multiprocessing import Pool
from tqdm import tqdm

from params import read_params_file
from properties import get_max_length
from metadata import get_metadata
import new_trench_algos

"""
Write two PyTables tables:
1. The trench row co-ordinates
2. The trench co-ordinates

Uses a peak detection algorithm to detect features (rows or trenches) after applying an unsharp mask to the image.
Typically, the idea is to use phase contrast images, in which these features stand out.

The row or trench co-ordinates can then be used as "regions" during segmentation analysis.

TODO: come up with a way to pass in which detection algorithm to use.
Add the ability to specify whether a given unit is in pixels or in microns in the YAML params file.
"""


def make_row_coords_type(max_filename_len):
    """
    Define the column types for the PyTable: this stores the coordinates of a trench row.
    Assumed to be shared across Z-levels.
    """
    # Define the PyTables column types using a dictionary
    column_types = {
        "info_file": tables.StringCol(max_filename_len),
        "info_fov": tables.UInt16Col(),
        "info_frame": tables.UInt16Col(),
        "info_z_level": tables.UInt16Col(),
        "info_row_number": tables.UInt16Col(),
        # min_row, min_col, max_row, max_col
        # Doesn't work with pandas
        # "bounding_box"    : tables.UInt16Col(shape=(4))
        "min_row": tables.UInt16Col(),
        "min_col": tables.UInt16Col(),
        "max_row": tables.UInt16Col(),
        "max_col": tables.UInt16Col(),
    }

    return type("RowCoords", (tables.IsDescription,), column_types)


def make_trench_coords_type(max_filename_len):
    """
    Define the column types for the PyTable: this stores the coordinates of a trench in a given trench row.
    Assumed to be shared across Z-levels.
    """
    # Define the PyTables column types using a dictionary
    column_types = {
        "info_file": tables.StringCol(max_filename_len),
        "info_fov": tables.UInt16Col(),
        "info_frame": tables.UInt16Col(),
        "info_z_level": tables.UInt16Col(),
        "info_row_number": tables.UInt16Col(),
        "info_trench_number": tables.UInt16Col(),
        # min_row, min_col, max_row, max_col
        # Doesn't work with pandas
        # "bounding_box"       : tables.UInt16Col(shape=(4,))
        "min_row": tables.UInt16Col(),
        "min_col": tables.UInt16Col(),
        "max_row": tables.UInt16Col(),
        "max_col": tables.UInt16Col(),
    }

    return type("TrenchCoords", (tables.IsDescription,), column_types)


def compile_trenches_to_table(file_names, metadata, file_to_frames, out_dir):
    """
    After the trench rows & trenches have been detected, compile the results to an
    HDF5 table.
    TODO: if some regions were shared across frames etc., then need to add duplicate
    entries for those shared regions.
    """
    max_filename_len = get_max_length(file_names)

    # Create the H5 file for writing the table to
    out_file_path = os.path.join(out_dir, "regions_tables.h5")
    out_h5file = tables.open_file(out_file_path, mode="w")

    # Create the trench rows table
    RowCoords = make_row_coords_type(max_filename_len)
    table_rowcoords = out_h5file.create_table(
        "/",
        "row_coords",
        RowCoords,
        filters=tables.Filters(complevel=1, complib="zlib"),
    )
    table_rowcoords_row = table_rowcoords.row

    # Create the trenches table
    TrenchCoords = make_trench_coords_type(max_filename_len)
    table_trenchcoords = out_h5file.create_table(
        "/",
        "trench_coords",
        TrenchCoords,
        filters=tables.Filters(complevel=1, complib="zlib"),
    )
    table_trenchcoords_row = table_trenchcoords.row

    # Iterate the same dimensions as before
    for f in file_names:
        for fov in metadata[f]["fields_of_view"]:
            for frame in file_to_frames[f]:
                for z in metadata[f]["z_levels"]:
                    # FIXME this doesn't take advantage of any linking system
                    path_string = "{}/FOV_{}/Frame_{}/Z_{}_regions.h5".format(
                        f, fov, frame, z
                    )
                    in_file_path = os.path.join(out_dir, path_string)
                    regions_file = tables.open_file(in_file_path, mode="r")

                    # Open the array containing 1 or more row coords
                    rows = regions_file.get_node("/row_coords").read()

                    # Write the row coords to a table
                    for i, row in enumerate(rows):
                        table_rowcoords_row["info_file"] = f
                        table_rowcoords_row["info_fov"] = fov
                        table_rowcoords_row["info_frame"] = frame
                        table_rowcoords_row["info_z_level"] = z
                        table_rowcoords_row["info_row_number"] = i

                        # NOTE wanted to write a single column with 4 bbox values,
                        # but that doesn't work with Pandas :(
                        # table_rowcoords_row["bounding_box"] = row
                        table_rowcoords_row["min_row"] = row[0]
                        table_rowcoords_row["min_col"] = row[1]
                        table_rowcoords_row["max_row"] = row[2]
                        table_rowcoords_row["max_col"] = row[3]

                        table_rowcoords_row.append()

                    # Open the array containing 1 or more row coords
                    trenches = regions_file.get_node("/trench_coords").read()

                    # Write the trench coords to a table
                    # First, go through each row in the vlarray
                    for i, row in enumerate(trenches):
                        # Each trench in each row
                        for j, trench in enumerate(row):
                            table_trenchcoords_row["info_file"] = f
                            table_trenchcoords_row["info_fov"] = fov
                            table_trenchcoords_row["info_frame"] = frame
                            table_trenchcoords_row["info_z_level"] = z
                            table_trenchcoords_row["info_row_number"] = i
                            table_trenchcoords_row["info_trench_number"] = j

                            # Doesn't work with pandas
                            # table_trenchcoords_row["bounding_box"] = trench
                            table_trenchcoords_row["min_row"] = trench[0]
                            table_trenchcoords_row["min_col"] = trench[1]
                            table_trenchcoords_row["max_row"] = trench[2]
                            table_trenchcoords_row["max_col"] = trench[3]

                            table_trenchcoords_row.append()

                    regions_file.close()

    # Flush the tables; done writing!
    table_rowcoords.flush()
    table_trenchcoords.flush()

    out_h5file.close()


def main_detection_function(out_dir, in_file, num_cpu, params_file):
    """
    Main detection function.
    """
    # Define a dict of possible algos: input a name, output the function.
    # Add new algorithms here!
    algo_dict = {"algo_1": new_trench_algos.launch}

    in_file = os.path.join(in_file, "data.h5")
    h5file = tables.open_file(in_file, mode="r")

    # Iterate all the files & get their metadata
    # Loop the nodes immediately under Images to get the file names
    file_names = [i._v_name for i in h5file.list_nodes("/Images")]

    # Generate dict. of FOVs to process, for each ND2 file
    file_to_frames = {}

    # Progress bar
    total = 0

    # Metadata from each ND2 file
    metadata = {}
    for f in file_names:
        # Get metadata for each file
        n = h5file.get_node("/Metadata/{}".format(f))()
        m = get_metadata(n)

        # Total # of frames to be processed, for progress bar
        total += len(m["fields_of_view"]) * len(m["frames"]) * len(m["z_levels"])

        # Each time frame has its own detection (do not share the regions across time frames).
        # NOTE is the sorting necessary, or do they always come out sorted from the conversion program?
        file_to_frames[f] = sorted(m["frames"])

        # store in a dictionary
        metadata[f] = m

    h5file.close()

    # Set up a progress bar
    pbar = tqdm(total=total, desc="Frame #")

    # Read the parameters
    params = read_params_file(params_file)

    # Determine which function to use for detecting trenches
    detection_func = algo_dict[params["algorithm"]]

    # Run in parallel. Implementation details are left to the algorithm.
    with Pool(processes=num_cpu) as pool:
        for f in file_names:
            # Make a copy of params, & edit it to convert from microns to pixels
            # Must do this for every "file," because each could have its own px_mu value.
            params_copy = params
            px_mu = metadata[f]["pixel_microns"]

            # Go through "parameters"
            for k1, param_dict in params_copy["parameters"].items():
                # For each set of parameters:
                for k2, parameter in param_dict.items():
                    # Convert?
                    # Replace the (value, dimension) dict at this level with just a value.
                    if parameter["dimension"] == "microns":
                        param_dict[k2] = int(parameter["value"] / px_mu)
                    else:
                        param_dict[k2] = parameter["value"]

            # Call the detection function
            # It's up to the function how to iterate FOVs, etc. and how to write
            # data to disk (although the format should be standard such that the
            # table creation at the end still works).
            # We provide a multi-processing pool, the HDF5 input file,
            # the frames within this file, all necessary metadata &
            # parameters, and a progress bar.
            detection_func(
                f,
                file_to_frames[f],
                pool,
                metadata[f],
                params_copy,
                in_file,
                out_dir,
                pbar,
            )

        pool.close()
        pool.join()

    # Done!

    # Done recording detected trench rows & trenches to HDF5-backed arrays.
    # Go through all of those arrays & compile the information into a table,
    # for easier searching. Store the table in the same HDF5 file.
    # All subsequent TrenchCoat steps assume that this table exists.
    # It's *much* easier to search for trench coordinates using a table!
    # NOTE this does raise the question of whether there's a way to just
    # create the table directly in the first place.
    # Could use the strategy of making small tables first, and then merging
    # them at the end.
    compile_trenches_to_table(file_names, metadata, file_to_frames, out_dir)
