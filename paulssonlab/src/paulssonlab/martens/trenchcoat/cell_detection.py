#!/usr/bin/env python3

import numpy
import time
import tables
import os

from scipy.signal import find_peaks
import argparse
from multiprocessing import Pool

# ### Find gaps between cells in trenches, and write bounding boxes to a table ###


def parse_prominence_values(file_path):
    """
    Input file path, output dict relating channel names to prominence values for peak detection
    """
    channel_to_prominence = {}
    with open(file_path, "r") as filehandle:
        for line in filehandle:
            line = line.strip()
            elements = line.split("\t")
            channel_name = elements[0]
            prominence = int(elements[1])

            channel_to_prominence[channel_name] = prominence

    return channel_to_prominence


def load_img_data(image_node, crop_top, crop_bottom, x_min, x_max, bg):
    """
    Input an image, crope given x/y bounds, and subtract a background value
    """
    data = image_node[crop_top:crop_bottom, x_min:x_max]
    data = bg_correction(data, bg)
    return data


def bg_correction(image, bg):
    """
    Basic background correction which doesn't allow underflow on pixels < bg
    """
    return (((image > bg) * image) + ((image <= bg) * bg)) - bg


def invert_intensity(image):
    """
    Inverts pixels within the range 0 _ max. Useful for detecting troughs (opposite of peaks)
    """
    return image.max() - image


def make_bbox_type():
    """
    For storing the bounding box information into the table
    """
    # Define the PyTables column types using a dictionary
    column_types = {
        "info_fov": tables.UInt16Col(),
        "info_frame_number": tables.UInt16Col(),
        "info_trench_number": tables.UInt16Col(),
        "info_channel_name": tables.StringCol(
            16
        ),  # TODO: use an enum column type? 16 chars is prob. enough to store a channel name for now.
        "info_bbox_number": tables.UInt16Col(),
        "min_row": tables.UInt16Col(),
        "min_col": tables.UInt16Col(),
        "max_row": tables.UInt16Col(),
        "max_col": tables.UInt16Col(),
    }

    return type("Bounding_Box", (tables.IsDescription,), column_types)


def run_bbox_detection(
    fov_file_path,
    identifier,
    crop_top,
    crop_bottom,
    y_dimension,
    prominence_file_path,
    min_peak_distance,
):
    """
    For a given FOV, scan each segmentation channel and search for dips in intensity. These probably correspond to gaps
    between cells. Write bounding boxes to a table which delineate the regions between gaps.
    """
    # Load prominence values
    prominence = parse_prominence_values(prominence_file_path)

    # Set seg. channels
    seg_channels = prominence.keys()

    # TODO: pass in this value as an argument; & maybe there are better ways to handle background removal?
    background = 100

    # Open this FOV
    fov_h5file = tables.open_file(fov_file_path, mode="r+")

    # Extract the FOV number from the file name
    base_name = os.path.basename(fov_file_path)
    (root, extension) = os.path.splitext(base_name)
    fov = int(root.split("_")[1])

    # Create new table for storing bounding box information, across all trench numbers & frames, in this FOV (will later merge these tables)
    Bounding_Box = make_bbox_type()

    table_name = "inner_trench_bounding_boxes_{}".format(identifier)

    # Create the bbox table
    # Create bbox tables for each FOV, then merge them together at the end
    if fov_h5file.__contains__("/tables/{}".format(table_name)):
        node = fov_h5file.get_node("/tables/{}".format(table_name))
        node._f_remove()

    bbox_table = fov_h5file.create_table(
        "/tables", table_name, Bounding_Box, "Inner_trench_bounding_boxes"
    )
    bbox_row = bbox_table.row

    # Input trench coordinates
    # Open the detected trenches table
    trenches_table_name = "/tables/trench_coordinates_{}".format(identifier)
    trenches_table = fov_h5file.get_node(trenches_table_name)

    ranges = [
        (row["bounding_box_min_col"], row["bounding_box_max_col"])
        for row in trenches_table.read()
    ]

    # For each frame
    for time_node in fov_h5file.root.images._f_iter_nodes():
        frame_number = int(time_node._v_name.split("_")[1])
        print("\t{}".format(time_node._v_name))

        # NOTE: would it be faster to first low the whole image within [crop_top : crop_bottom],
        # and then afterwards further crop this to get each trench? (Fewer calls back to the HDF5 file, vs. more total data to transfer)
        # TODO: not specify max
        # img_row = load_img_data(time_node._f_get_child(c), crop_top, crop_bottom, 0, 2048, background)

        # Loop all the xmin, xmax pairs
        for n, (xmin, xmax) in enumerate(ranges):

            # Input channel, output trench image data
            for c in seg_channels:

                # Load the data
                img = load_img_data(
                    time_node._f_get_child(c),
                    crop_top,
                    crop_bottom,
                    xmin,
                    xmax,
                    background,
                )
                # img = img_row[:, xmin : xmax]

                # Invert the intensities, so that we can detect space between cells as peaks
                inv = invert_intensity(img)

                # Calculate the mean in the x-dimension
                avg = inv.mean(axis=1)

                # Run peak detection in the y-dimension
                peaks, _ = find_peaks(
                    avg, prominence=prominence[c], distance=min_peak_distance
                )

                # No peaks were detected, so the box is number zero and spans entire y range
                if len(peaks) == 0:
                    bbox_row["info_fov"] = fov
                    bbox_row["info_frame_number"] = frame_number
                    bbox_row["info_trench_number"] = n
                    bbox_row["min_row"] = xmin
                    bbox_row["max_row"] = xmax
                    bbox_row["info_channel_name"] = c
                    bbox_row["info_bbox_number"] = 0
                    bbox_row["min_col"] = 0
                    bbox_row["max_col"] = y_dimension

                    # Append bounding box information to the table
                    bbox_row.append()

                # Iterate each detected peak
                else:
                    for i, p in enumerate(peaks):
                        bbox_row["info_fov"] = fov
                        bbox_row["info_frame_number"] = frame_number
                        bbox_row["info_trench_number"] = n
                        bbox_row["min_row"] = xmin
                        bbox_row["max_row"] = xmax
                        bbox_row["info_channel_name"] = c
                        bbox_row["info_bbox_number"] = i

                        if i == 0:
                            bbox_row["min_col"] = 0
                            bbox_row["max_col"] = peaks[i]

                        elif i > 0:
                            bbox_row["min_col"] = peaks[i - 1]
                            bbox_row["max_col"] = peaks[i]

                        # Append bounding box information to the table
                        bbox_row.append()

                    # And also write the remaining space after the final peak
                    # TODO: what if this space is empty? Do we care?
                    bbox_row["info_fov"] = fov
                    bbox_row["info_frame_number"] = frame_number
                    bbox_row["info_trench_number"] = n
                    bbox_row["min_row"] = xmin
                    bbox_row["max_row"] = xmax
                    bbox_row["info_channel_name"] = c
                    bbox_row["info_bbox_number"] = i
                    bbox_row["min_col"] = peaks[i]
                    bbox_row["max_col"] = y_dimension

                    # Append bounding box information to the table
                    bbox_row.append()

    # Done!
    bbox_table.flush()


# TODO: throw errors if the arguments don't make sense?
def parse_args():
    parser = argparse.ArgumentParser(
        description="Find gaps between cells in trenches, and write bounding boxes to a table."
    )
    parser.add_argument("-i", "--in-dir", type=str, help="Input HDF5 directory path")
    parser.add_argument(
        "-p", "--parent", type=str, help="Parent HDF5 file (with links & merged tables)"
    )
    parser.add_argument("-n", "--num-cpu", type=int, help="Number of CPUs to use")
    parser.add_argument("-I", "--identifier", type=str, help="Trench row identifier")
    parser.add_argument("-T", "--crop-top", type=int, help="Crop top")
    parser.add_argument("-B", "--crop-bottom", type=int, help="Crop bottom")
    parser.add_argument(
        "-d", "--min-peak-distance", type=int, help="Minimum peak distance"
    )
    parser.add_argument(
        "-P",
        "--prominence-file",
        type=str,
        help="File with prominence values for each channel",
    )

    return parser.parse_args()


###

if __name__ == "__main__":
    ### Parse command-line arguments
    args = parse_args()

    in_dir = args.in_dir
    parent_file = args.parent
    num_cpu = args.num_cpu
    identifier = args.identifier
    crop_top = args.crop_top
    crop_bottom = args.crop_bottom
    min_peak_distance = args.min_peak_distance
    prominence_file_path = args.prominence_file

    # FIXME can't send in a dict through the starmap?
    # -> Have to re-read & parse these values every single time... there must be a better way to do this.
    # One way could be parallel lists or something, and then convert to a dict internally.
    ## Load prominence values
    # prominence = parse_prominence_values(prominence_file_path)

    ## Set seg. channels
    # seg_channels = prominence.keys()

    # FIXME what's the right way to read these in?
    y_dimension = crop_bottom - crop_top

    ### Find the bounding boxes & write to tables for each FOV
    print("Detecting cell boundaries...")
    start = time.time()

    fov_dir = os.path.join(in_dir, "FOV")
    files = os.listdir(fov_dir)
    ## DEBUG
    # files = ["FOV_0.h5"]

    # Run in parallel
    if num_cpu > 1:
        args = [
            (
                os.path.join(fov_dir, filename),
                identifier,
                crop_top,
                crop_bottom,
                y_dimension,
                prominence_file_path,
                min_peak_distance,
            )
            for filename in files
        ]

        # TODO: weed out any files which do not end in h5

        with Pool(processes=num_cpu) as p:
            p.starmap(run_bbox_detection, args)

    # Do not run in parallel
    else:
        for filename in files:
            if filename.endswith("h5"):
                run_bbox_detection(
                    os.path.join(fov_dir, filename),
                    identifier,
                    crop_top,
                    crop_bottom,
                    y_dimension,
                    prominence_file_path,
                    min_peak_distance,
                )

    ### Done looping all FOV
    end = time.time()
    print(end - start)
    print("Done detecting cell boundaries.")

    ### Merge the tables
    print("Merging tables...")
    start = time.time()

    BBox = make_bbox_type()

    h5file = tables.open_file(os.path.join(in_dir, parent_file), mode="r+")
    table_name = "inner_trench_bounding_boxes_{}".format(identifier)

    # Create a table for storing the properties of the entire trench
    # If table exists, delete it and start over
    if h5file.__contains__("/tables/{}".format(table_name)):
        node = h5file.get_node("/tables/{}".format(table_name))
        node._f_remove()

    bbox_table_global = h5file.create_table(
        h5file.root.tables, table_name, BBox, "Inner_trench_bounding_boxes"
    )

    for filename in os.listdir(fov_dir):
        if filename.endswith("h5"):
            table_file = tables.open_file(os.path.join(fov_dir, filename), "r")
            table = table_file.get_node("/tables/{}".format(table_name))

            for next_row in table.iterrows():
                bbox_table_global.append([next_row[:]])

            table.close()
            table_file.close()

    # Done copying trench properties
    bbox_table_global.flush()
    bbox_table_global.close()
    h5file.close()

    end = time.time()
    print(end - start)
    print("Done merging tables.")

    ### Done!
