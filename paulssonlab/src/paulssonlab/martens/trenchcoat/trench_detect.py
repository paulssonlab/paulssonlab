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
from scipy.signal import savgol_filter

"""
Write two PyTables tables:
1. The trench row co-ordinates
2. The trench co-ordinates

Uses a peak detection algorithm to detect features (rows or trenches) after applying an unsharp mask to the image.
Typically, the idea is to use phase contrast images, in which these features stand out.

The row or trench co-ordinates can then be used as "regions" during segmentation analysis.
"""


def find_trenches(
    img,
    axis=1,
    crop_top=0,
    crop_bottom=0,
    unsharp_mask_radius=5,
    unsharp_mask_amount=30.0,
    cutoff=0.10,
    plateau_size=4,
    distance=15,
    padding_ratio=0.5,
    window_length=5,
    polyorder=1,
):
    """
    Find trenches within an image (or image slice). Ideally, the image slice only contains a single row of trenches.

    Cropping can help, if the trench row isn't perfectly aligned,
    by removing regions outside of the trenches area.

    The unsharp mask filter helps to emphasize the edges of trenches.

    The cutoff is used to separate trenches from inter-trench space.

    The distance helps filter out the non-trenches, since their peaks are less than those of trenches.

    The padding ratio helps to add back the darker edges of the trenches, which don't make the cutoff.
    """
    # The unsharp mask helps to accentuate local differences by "sharpening" the image.
    # NOTE unsharp mask gives worse results when operating on floats. Not sure why?
    unsharp = unsharp_mask(
        img.astype(numpy.uint16), radius=unsharp_mask_radius, amount=unsharp_mask_amount
    )

    # The mean is useful as a normalization
    mean = numpy.mean(unsharp, axis=axis)

    # Run smoothing: this helps merge twin peaks from the edges
    # of trenches with cells, and also the twin peaks for the space
    # in between trenches.
    sg = savgol_filter(mean, window_length=window_length, polyorder=polyorder)

    # Define an arbitrary threshold, below which we ignore,
    # and above which we deem peak-worthy.
    gt = sg > cutoff

    # By defining the plateau size as such, we can eliminate the peaks which arise
    # from the emblazoned numerals. Plateau size helps to eliminate very small
    # detected peaks, which are actually just noise.
    #
    # Distance is crucial in filtering out the "peaks" arising from the inter-trench
    # spacing. These inter-trench spaces are almost always dimmer than trenches,
    # even with cells in them.
    #
    # Note that trenches with cells generate two peaks, with a dip in the middle
    # where the cells are, but these two peaks are merged as long as the entirety
    # (peak, drop, peak) is greater than the low cutoff.
    peaks, peaks_properties = find_peaks(
        gt, plateau_size=plateau_size, distance=distance
    )

    # If it fails to find any peaks, then end now & return an empty array
    if peaks.size == 0:
        print("Warning: no peaks found!")
        return peaks

    # Then, it will be important to pad the left_edges and right_edges a bit,
    # since the trench edges are dark and don't meet the cutoff, but can
    # contain a lot of fluorescence intensity. To do this, define a padding,
    # and subtract it from all left_edges, and add it to all right edges.
    median_detected_width = int(
        numpy.median(peaks_properties["right_edges"] - peaks_properties["left_edges"])
    )

    # Apply padding as some fraction or multiple of the median trench width.
    # So first determine a uniform width, perhaps the median width,
    # and then add +/- 50% to left & right sides.
    # NOTE: With a value of 1.0, this will triple the effective trench width.
    padding = int(padding_ratio * median_detected_width)

    # Adjust all trenches to have equal widths
    # NOTE: is the max value img.shape[0] or [1]?
    (left_edges, right_edges) = correct_widths(
        peaks_properties["left_edges"],
        peaks_properties["right_edges"],
        median_detected_width,
        padding,
        img.shape[0],
    )

    # Return rectangular coordinates for the left and right edges, sharing the min and max col.
    result = numpy.array(
        [[i[0], 0, i[1], img.shape[1]] for i in zip(left_edges, right_edges)]
    )

    return result


def correct_widths(
    left_edges_old, right_edges_old, median_detected_width, padding, max_value
):
    """
    If a given trench is not equal to the median dimension, then it is either too wide, or too narrow.
    The difference between the trench and the median is either odd, or even.
    If the difference is even, then add (or subtract) half of this difference from the left & the right.
    If the difference is odd, then we have to decide whether to add/remove 1 extra pixel from the left,
    or from the right.
    Let's just (arbitrarily decide to give more to the right.
    """
    left_edges = []
    right_edges = []

    # FIXME the width correction procedure doesn't work!
    for l, r in zip(left_edges_old, right_edges_old):
        difference = median_detected_width - (r - l)
        left = l - padding
        right = r + padding

        # Wider than the median, must remove some width
        if difference < 0:
            if difference % 2 == 0:
                left -= difference // 2
                right += difference // 2
            else:
                # Integer divide with negative numbers!
                # https://stackoverflow.com/questions/19517868/integer-division-by-negative-number
                left -= ((-1 * difference) // 2) * -1
                right += ((-1 * difference) // 2) * -1 - 1

        # Narrower than the median, must add some width
        elif difference > 0:
            if difference % 2 == 0:
                left -= difference // 2
                right += difference // 2
            else:
                left -= difference // 2
                right += difference // 2 + 1

        # NOTE that now it is possible for a trench to extend outside of the image.
        # Filter out any peaks which have left_edge < 0 or right_edge > img_dimension.
        # by not adding them to the final list.
        if (left > 0) & (right < max_value):
            left_edges.append(left)
            right_edges.append(right)

    # All trenches now have padding and uniform widths and remain within the image bounds.

    return (left_edges, right_edges)


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

    # DEBUG
    # img = numpy.rot90(img, -1)
    # params["trench_rows"]["axis"] = 0
    # params["trench_rows"]["axis_two"] = 1

    # 2. Create an HDF5 file for writing trench row & trench regions info.
    dir_path = os.path.join(out_dir, "{}/FOV_{}/Frame_{}/".format(filename, fov, frame))
    pathlib.Path(dir_path).mkdir(parents=True, exist_ok=True)

    out_file_path = os.path.join(dir_path, "Z_{}_regions.h5".format(z_level))
    regions_file = tables.open_file(out_file_path, mode="w")

    # 3. Analyze the image to find trench rows & trenches, write coordinates of the detected rows & trenches to HDF5
    rows = detect_trench_rows(regions_file, params, img, filename, fov, frame, z_level)
    print(rows)
    # DEBUG
    # print(rows.shape)

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


def detect_trench_rows(regions_file, img, axis, cutoff, plateau_size):
    """
    Analyze the image to find trench rows, write coordinates to HDF5
    """
    # Peak detection parameters for rows of trenches: find all the trench rows in the image.
    if "trench_rows" in params:
        # The variance is expected to be greater where there are features (numerals or trenches)
        var = numpy.var(img, axis=axis)

        # The mean is useful as a normalization
        mean = numpy.mean(img, axis=axis)

        # Normalize the variance by the mean
        var_mean = var / mean

        # Define an arbitrary threshold, below which we ignore,
        # and above which we deem peak-worthy
        gt = var_mean > cutoff

        # By defining the plateau size as such, we can eliminate the peaks
        # which arise from the emblazoned numerals.
        peaks, peaks_properties = find_peaks(gt, plateau_size=plateau_size)

        # FIXME img.shape[0] or [1] ? (width, or height?)
        rows = numpy.array(
            [
                [0, i[0], img.shape[1], i[1]]
                for i in zip(
                    peaks_properties["left_edges"], peaks_properties["right_edges"]
                )
            ]
        )

        # Adjust the coords to the raw image (pre-cropped)
        for r in rows:
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
        trenches = find_trenches(
            img[row[0] : row[2], row[1] : row[3]], **params["trenches"]
        )

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
    # FIXME this function does both rows & trenches ~ but what if we only have rows, but no trenches?
    # Should we even allow row detection without trench detection? What would be the point?
    # if "trenches" in params_copy:
    # compile_trenches_to_table(file_names, metadata, file_to_frames, out_dir)
