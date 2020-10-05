#!/usr/bin/env python3

import os
from tqdm import tqdm
import numpy
from scipy.signal import find_peaks
from skimage.filters import unsharp_mask
import tables
from multiprocessing import Pool
from params import read_params_file
from properties import get_max_length
from metadata import get_metadata
import pathlib

"""
Write two PyTables tables:
1. The trench row co-ordinates
2. The trench co-ordinates

Uses a peak detection algorithm to detect features (rows or trenches) after applying an unsharp mask to the image.
Typically, the idea is to use phase contrast images, in which these features stand out.

The row or trench co-ordinates can then be used as "regions" during segmentation analysis.
"""

# TODO what needs to change if the images are now stored with F-ordering?


def run_trench_analysis(in_file, out_dir, filename, fov, frame, z_level, params):
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
    rows = detect_trench_rows(regions_file, params, img, filename, fov, frame, z_level)

    # 4. For each row, find all the trenches
    if "trenches" in params:
        # FIXME this is still the cropped image, but the co-ordinates are for uncropped ~ 2 solutions:
        # a. add back the crop dimensions,
        # b. re-load the image w/o cropping it
        detect_trenches_in_rows(
            filename, fov, frame, z_level, rows, img, params, regions_file
        )

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
            print("Error loading cropped image. Invalid node? Missing crop params?")
            pass
    else:
        # Load the whole image
        # TODO is there a way to get it to read into an empty numpy array with F-ordering?
        # Otherwise, we either have to rotate the image, or flip all of the co-ordinates.
        img = h5file.get_node(node_path).read()

        # And define the missing crop params, so that they are always defined (for convenience)
        # FIXME width & height?
        params["crop"] = {"top": 0, "bottom": img.height, "left": 0, "right": img.width}

    h5file.close()

    return img


def detect_trench_rows(regions_file, params, img, filename, fov, frame, z_level):
    """
    Analyze the image to find trench rows & trenches, write coordinates of the detected rows & trenches to HDF5
    """
    # Peak detection parameters for rows of trenches: find all the trench rows in the image.
    if "trench_rows" in params:
        rows = detect_rows_or_trenches(img, params["trench_rows"])  # , params["debug"])

        # Adjust the coords to the raw image (pre-cropped)
        for r in rows:
            # DEBUG
            # print(r)
            r[0] -= params["crop"]["left"]
            r[1] -= params["crop"]["top"]
            r[2] -= params["crop"]["left"]
            r[3] -= params["crop"]["top"]

    # If no trench row detection parameters were specified, then there is just a single row,
    # equal to the original cropping parameters.
    # If no cropping is desired, then the params must span the whole image.
    # Return a 2D array with a single entry, so that the dims are compatible
    # for when either just 1 row is detected, or more than 1 are detected.
    else:
        rows = numpy.array(
            [
                [
                    params["crop"]["left"],
                    params["crop"]["top"],
                    params["crop"]["right"],
                    params["crop"]["bottom"],
                ]
            ]
        )

    # Write all row coords to disk as a single array
    regions_file.create_array("/", "row_coords", obj=rows)

    return rows


def detect_trenches_in_rows(
    filename, fov, frame, z_level, rows, img, params, regions_file
):
    # Store all detected trenches as a list of numpy arrays, which will then be stored as a variable-length pytables array.
    # NOTE: reading from the vlarray will return a list of numpy arrays, each with its own length.
    detected_trenches_list = []

    for i, row in enumerate(rows):
        # Detect trenches within the row
        regions = detect_rows_or_trenches(img, params["trenches"])  # , params["debug"])

        # Convert coordinates back to pre-cropped (not in the region, but in the raw image).
        # NOTE confused. Why is it 0->1, 1->0, 2->1, 3->0 and not the converse?
        for r in regions:
            r[0] += row[1]
            r[1] += row[0]
            r[2] += row[1]
            r[3] += row[0]

        # Append to the global list
        detected_trenches_list.append(regions)

    # Write out the trench coordinates using a variable-length array in HDF5
    results = regions_file.create_vlarray(
        "/",
        "trench_coords",
        atom=tables.UInt16Atom(shape=(4)),
        expectedrows=len(detected_trenches_list),
    )

    for i in detected_trenches_list:
        results.append(i)

    # Done!
    results.close()


def detect_rows_or_trenches(img, params):
    """
    Find rows of trenches, or trenches within a row.

    Input a single image (which may have already been cropped), and parameters:
    1. axis along which to do the flattening (x for detecting rows, y for detecting trenches)
    2. unsharp_mask filter: radius, amount
    3. min (& max?) cutoffs
    4. min distance between features

    Return a 3D numpy array:
    Dim 1: # of the feature (row0,row1,... or trench0,trench1,...)
    Dim 2: [min_row, min_col, max_row, max_col] values in the *cropped* image.
    [requires then editing to get raw coordinates.]

    , each row containing 4 co-ordinates:
    min_row, min_col, max_row, max_col, defining the rectangular region corresponding to that row of trenches.

    TODO sanity checking on the result, to make sure that it's within the boundaries of the image?
    FIXME what are the necessary parameters?
    """
    # Run the unsharp mask to enhance edges (e.g. trenches in phase contrast) & to normalize from 0.0 - 1.0
    # Values that worked for trenches: radius=2, amount=40
    img = unsharp_mask(
        img, radius=params["unsharp_mask_radius"], amount=params["unsharp_mask_amount"]
    )

    # Flatten in the x-dimension: 1, y-dimension: 0
    peaks_data = img.mean(axis=params["axis"])

    # Dimension of the opposite axis (e.g. Width of image when detecting trenches)
    img_dim = img.shape[params["axis_two"]]

    # Use scipy's peak detection algorithm
    peaks, _ = find_peaks(peaks_data)

    # Remove non-trench peaks
    feature_peaks = remove_non_features(peaks_data, peaks, params["min_cutoff"])

    # Merge peaks which are close to each other
    merged_peaks = merge_peaks(feature_peaks, params["min_distance"])

    # Starting, ending x-coordinates for each trench.
    # These can be used to define rectangular regions of interest.
    # FIXME converting trench width to half width? best way?
    # If it's always the same, then pass this in rather than the whole width?
    feature_half_dim = params["feature_dimension"] // 2
    # FIXME the max_col param should be the largest possible indexing value for trench heights
    # -> equal to the cropped upper bounds.
    ranges = define_features(merged_peaks, feature_half_dim, 0, img.shape[1])

    # TODO 2 cleanup steps:

    # 1. Filter out overlapping features
    # UNIMPLEMENTED!

    # 2. Remove out-of-bounds features
    # TODO: replace if statements with while() loops, in case there are multiple trenches to be removed?
    # (if they are overlapping...)

    # Remove first or last trench if they are out of bounds
    if len(ranges) != 0:
        if ranges[0][0] < 0:
            ranges.pop(0)

    if len(ranges) != 0:
        if ranges[-1][2] > img_dim:
            ranges.pop()

    # Convert to numpy array
    return numpy.array(ranges)


def remove_non_features(peaks_data, peaks, min_cutoff):
    """
    Remove non-features using min_cutoff.
    """
    # Filter out: only keep inter-feature peak indices
    threshold = (peaks_data[peaks] > min_cutoff) * peaks

    # Remove zero indices
    nz = numpy.nonzero(threshold)

    return peaks[nz]


def merge_peaks(peaks, min_peak_distance):
    """
    Merge the nearby peaks in peaks, a numpy array resulting from
    find_peaks().

    TODO: also add a sanity check: only merge them if both are above or below the threshold line
    """
    merged_peaks = []

    # Use double while loops to allow for merging clusters of nearby peaks
    i = 0
    while i < len(peaks) - 1:
        cluster = []
        cluster.append(peaks[i])

        j = i + 1
        # NOTE: comparison is always made to the i'th peak. This might behave weirdly
        # if there are many peaks, in series, which are spaced ~= min_peak_distance from one another
        while (j < len(peaks)) and (peaks[j] - peaks[i] < min_peak_distance):
            cluster.append(peaks[j])
            j += 1

        avg = numpy.mean(cluster).astype(numpy.int)
        merged_peaks.append(avg)
        i = j

    # Finally, the last peak
    # If i > len(peaks) - 1, then we merged the last one, and so can skip it,
    # but if i == len(peaks) - 1, then we have to keep it.
    if i == len(peaks) - 1:
        merged_peaks.append(peaks[i])

    # Convert from numpy list to numpy array
    return numpy.array(merged_peaks)


def define_features(peaks, feature_half_dimension, min_col, max_col):
    """
    Detect Features using just the feature midpoints & the pre-specified dimensions:
    define left -> right boundaries for each trench, or
    top -> bottom boundaries for each row.

    FIXME Assume each trench spans from top to bottom. --> change to not assuming the other dims.

    TODO: sanity bounds check near the edges of the FOV (when adding or subbing trench_half_width)
    Before making ranges, check whether would overlap?

    TODO: return 4 co-ordinates: either
    min_row, min_col, max_row, max_col,
    OR x_min, x_max, y_min, y_max
    """
    ranges = [
        (t - feature_half_dimension, min_col, t + feature_half_dimension, max_col)
        for t in peaks
    ]

    return ranges


def compile_trenches_to_table(file_names, metadata, file_to_frames, out_dir):
    max_filename_len = get_max_length(file_names)

    # Create the H5 file for writing the table to
    out_file_path = os.path.join(out_dir, "regions_tables.h5")
    out_h5file = tables.open_file(out_file_path, mode="w")

    # Create the trench rows table
    RowCoords = make_row_coords_type(max_filename_len)
    table_rowcoords = out_h5file.create_table("/", "row_coords", RowCoords)
    table_rowcoords_row = table_rowcoords.row

    # Create the trenches table
    TrenchCoords = make_trench_coords_type(max_filename_len)
    table_trenchcoords = out_h5file.create_table("/", "trench_coords", TrenchCoords)
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

                        # Doesn't work with pandas
                        # table_rowcoords_row["bounding_box"]    = row
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
                            # table_trenchcoords_row["bounding_box"]       = trench
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

            if "trenches" in params_copy:
                params_copy["trenches"]["feature_dimension"] = int(
                    params_copy["trenches"]["feature_dimension"] / px_mu
                )
                params_copy["trenches"]["min_distance"] = int(
                    params_copy["trenches"]["min_distance"] / px_mu
                )

            if "trench_rows" in params_copy:
                params_copy["trench_rows"]["feature_dimension"] = int(
                    params_copy["trench_rows"]["feature_dimension"] / px_mu
                )
                params_copy["trench_rows"]["min_distance"] = int(
                    params_copy["trench_rows"]["min_distance"] / px_mu
                )

            for fov in metadata[f]["fields_of_view"]:
                for frame in file_to_frames[f]:
                    for z in metadata[f]["z_levels"]:
                        args = [in_file, out_dir, f, fov, frame, z, params_copy]
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
    compile_trenches_to_table(file_names, metadata, file_to_frames, out_dir)
