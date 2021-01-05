#!/usr/bin/env python3

import os
import pathlib

import numpy
import tables
from tqdm import tqdm
from multiprocessing import Pool

from params import read_params_file
from properties import get_max_length
from metadata import get_metadata
from new_trench_algos import find_trenches, generate_boxes

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


def run_trench_analysis(in_file, out_dir, filename, fov, frame, z_level, params, boxes):
    """
    Run basic trench analysis on an HDF5 file: detect trenches for each time
    frame, or just using the first time frame. Then, measure trench properties
    in all time frames.

    Store detected trenches in a PyTables array.

    Each array entry (row) has 4 values (min_row, min_col, max_row, max_col), and there is 1 row per region (trench).
    Regions must be rectangular, and must have the same dimension (for stacking purposes).
    """
    # 1. Load the image, optionally crop it
    img = load_image_for_regions_analysis(
        in_file, filename, fov, frame, z_level, params
    )

    # 2. Create an HDF5 file for writing trench row & trench regions info.
    dir_path = os.path.join(out_dir, "{}/FOV_{}/Frame_{}/".format(filename, fov, frame))
    pathlib.Path(dir_path).mkdir(parents=True, exist_ok=True)

    out_file_path = os.path.join(dir_path, "Z_{}_regions.h5".format(z_level))
    regions_file = tables.open_file(out_file_path, mode="w")

    # 3. Analyze the image to find trench rows & trenches, write coordinates of the detected rows & trenches to HDF5

    ## 4. For each row, find all the trenches
    # if "trenches" in params:
    ## FIXME this is still the cropped image, but the co-ordinates are for uncropped ~ 2 solutions:
    ## a. add back the crop dimensions,
    ## b. re-load the image w/o cropping it
    # detect_trenches_in_rows(filename, fov, frame, z_level, rows, img, params, regions_file)

    # UPDATE NEW trench row & trench detection algorithm
    # This detects both the rows & trenches in a single step
    # FIXME no correction for cropping; and the results are NOT written to disk!
    # To do so, need to retroactively parse out the row coords from the trenches.

    (trenches, y_offsets) = find_trenches(img, boxes, **params["trenches"])

    # Take the y_offsets and create bounding boxes for the rows of trenches
    rows = []
    for y in y_offsets:
        bbox = [0, y, img.shape[0], y + params["trenches"]["tr_height"]]
        rows.append(bbox)

    # Write all row coords to disk as a single array
    regions_file.create_array("/", "row_coords", obj=rows)

    # Write out the trench coordinates using a variable-length array in HDF5
    results = regions_file.create_vlarray(
        "/",
        "trench_coords",
        atom=tables.UInt16Atom(shape=(4)),
        expectedrows=len(trenches),
    )

    # Iterate the rows
    for r in trenches:
        results.append(r)

    # Done!
    regions_file.close()


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


def load_image_for_regions_analysis(in_file, filename, fov, frame, z_level, params):
    """
    Open the image to be used for features detection, & optionally apply cropping.
    """
    h5file = tables.open_file(in_file, mode="r")
    node_path = "/Images/{}/FOV_{}/Frame_{}/Z_{}/{}".format(
        filename, fov, frame, z_level, params["channel"]
    )

    # NOTE image was written in F-order. When reading from disk, image is "rotated" 90 degrees
    # counter-clockwise. For best performance, it makes sense to not re-rotate the image,
    # but instead to accept this new reference frame & interpret left/right/top/bottom,
    # and axis 0,1 appropriately.

    if "crop" in params:
        try:
            img = h5file.get_node(node_path)[
                params["crop"]["top"] : params["crop"]["bottom"],
                params["crop"]["left"] : params["crop"]["right"],
            ]
        except:
            # TODO Error!
            print(
                "Error loading cropped image. Invalid node? Missing crop params? Check channel name?"
            )
            pass
    else:
        # Load the whole image
        # TODO is there a way to get it to read into an empty numpy array with F-ordering?
        # Otherwise, we either have to rotate the image, or flip all of the co-ordinates.
        img = h5file.get_node(node_path).read()

        # And define the missing crop params, so that they are always defined (for convenience)
        # FIXME width & height? 0 & 1, or 1 & 0?
        params["crop"] = {
            "top": 0,
            "bottom": img.shape[0],
            "left": 0,
            "right": img.shape[1],
        }

    h5file.close()

    return img


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
        "/", "trench_coords", TrenchCoords, tables.Filters(complevel=1, complib="zlib")
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


def main_detection_function(out_dir, in_file, num_cpu, params_file, share_regions):
    """
    Main detection function.
    """
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
    # TODO re-work how region sharing would be handled:
    # Must specify which dims. to share across,
    # and then must duplicate the table entries at the end!
    if share_regions:
        for f in file_names:
            # Get metadata for each file
            n = h5file.get_node("/Metadata/{}".format(f))()
            m = get_metadata(n)

            # Total # of frames to be processed, for progress bar
            total += len(m["fields_of_view"]) * len(m["z_levels"])

            # Use the lowest numbered time frame to detect trenches in all time frames
            # (share the region across time frames).
            # NOTE is the sorting necessary, or do they always come out sorted from the conversion program?
            # Make it a list, so that the code in the main loop still works
            # file_to_frames[f] = [sorted(m["frames"][0])]
            # FIXME!! Temporarily just set it to zero.
            # TODO allow the user to specify which frame to use
            file_to_frames[f] = [0]

            # store in a dictionary
            metadata[f] = m
    else:
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

    # Read the parameters, and then convert the distances using pixel microns
    params = read_params_file(params_file)

    pbar = tqdm(total=total, desc="Frame #")

    def update_pbar(*a):
        pbar.update()

    # Run in parallel. Each File, FOV or Frame, can be processed independently.
    with Pool(processes=num_cpu) as p:
        for f in file_names:
            # Make a copy of params, & edit it to convert from microns to pixels
            params_copy = params
            px_mu = metadata[f]["pixel_microns"]

            # FIXME need a more generic way to convert specific params to pixels from microns
            # How to flag them for conversion?
            if "trenches" in params_copy:
                # params_copy["trenches"]["feature_dimension"] = int(params_copy["trenches"]["feature_dimension"] / px_mu)
                # params_copy["trenches"]["min_distance"]      = int(params_copy["trenches"]["min_distance"]      / px_mu)
                params_copy["trenches"]["tr_height"] = int(
                    params_copy["trenches"]["tr_height"] / px_mu
                )

            if "trench_rows" in params_copy:
                params_copy["trench_rows"]["feature_dimension"] = int(
                    params_copy["trench_rows"]["feature_dimension"] / px_mu
                )
                params_copy["trench_rows"]["min_distance"] = int(
                    params_copy["trench_rows"]["min_distance"] / px_mu
                )

            # Generate a set of rectangles which represent trenches
            # but which do not yet have known x,y starting points
            # which best match one or more of the trench rows.
            # Re-use the same boxes for all iterations!
            # TODO: more generic way to pass in shared variables,
            # in an algorithm-agnostic way?
            # Note that performance seems the same with or without
            # re-computing the boxes. But if they really are the same,
            # it makes sense to only calculate them once.
            # Stick boxes into params_copy somehow?
            boxes = generate_boxes(
                params_copy["trenches"]["tr_width"],
                params_copy["trenches"]["tr_height"],
                params_copy["trenches"]["tr_spacing"],
                metadata[f]["width"],
            )

            for fov in metadata[f]["fields_of_view"]:
                for frame in file_to_frames[f]:
                    for z in metadata[f]["z_levels"]:
                        args = [in_file, out_dir, f, fov, frame, z, params_copy, boxes]
                        p.apply_async(run_trench_analysis, args, callback=update_pbar)
                        # run_trench_analysis(*args) # DEBUG

        p.close()
        p.join()

    # Done looping all files
    pbar.close()

    # Done recording detected trench rows & trenches to HDF5-backed arrays.

    # TODO Link all of the individual h5 files into a single file ?

    # Go through all of those arrays & compile the information into a table,
    # for easier searching. Store the table in the same HDF5 file.
    # FIXME this function does both rows & trenches ~ but what if we only have rows, but no trenches?
    # Should we even allow row detection without trench detection? What would be the point?
    if "trenches" in params_copy:
        compile_trenches_to_table(file_names, metadata, file_to_frames, out_dir)
