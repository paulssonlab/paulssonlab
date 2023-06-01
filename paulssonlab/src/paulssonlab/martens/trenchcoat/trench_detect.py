#!/usr/bin/env python3
import os
from multiprocessing import Pool, get_context

import new_trench_algos
import numpy
import tables
from metadata import get_metadata
from params import read_params_file
from properties import get_max_length
from tqdm import tqdm

"""
Write two PyTables tables:
1. The trench row co-ordinates
2. The trench co-ordinates

Uses a peak detection algorithm to detect features (rows or trenches) after applying an unsharp mask to the image.
Typically, the idea is to use phase contrast images, in which these features stand out.

The row or trench co-ordinates can then be used as "regions" during segmentation analysis.

TODO: come up with a way to pass in which detection algorithm to use.
Add the ability to specify whether a given unit is in pixels or in microns in the YAML params file.

FIXME did I delete the dir creation step?
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


def main_detection_function(out_file, in_file, num_cpu, params_file):
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

    ###
    max_filename_len = get_max_length(file_names)

    # Create the H5 file for writing the table to
    out_h5file = tables.open_file(out_file, mode="w")

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
    ###

    # Set up a progress bar
    pbar = tqdm(total=total, desc="Frame #")

    # Read the parameters
    params = read_params_file(params_file)

    # Determine which function to use for detecting trenches
    detection_func = algo_dict[params["algorithm"]]

    # Run in parallel. Implementation details are left to the algorithm.
    # See: https://pythonspeed.com/articles/python-multiprocessing/
    # for why we use get_context("spawn")
    with get_context("spawn").Pool() as pool:
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
            # parameters, a progress bar, and two tables to write results to.
            results = detection_func(
                f,
                file_to_frames[f],
                pool,
                metadata[f],
                params_copy,
                in_file,
                table_rowcoords_row,
                table_trenchcoords_row,
                pbar,
            )

        # Close the pool
        pool.close()
        pool.join()

        # Flush the tables; done writing!
        table_rowcoords.flush()
        table_trenchcoords.flush()
        out_h5file.close()

    # Done!
