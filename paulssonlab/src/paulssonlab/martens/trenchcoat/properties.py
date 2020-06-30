#!/usr/bin/env python3

import numpy
import tables

"""
Functions for handling ND2 cell properties (measurements), including HDF5 tables.
"""


def write_properties_to_table(
    name,
    fov,
    frame,
    z_level,
    region_set_number,
    region_number,
    properties,
    sc,
    params,
    row,
    stack,
    ch_to_index,
):
    """
    Once a given trench has a segmentation mask, then write the properties of the masked region to the table.
    NOTE future versions of skimage will allow spitting out a properties object all at once, rather than lazily calculating them one at a time.
    """
    # NOTE skimage is changing their coordinates -- do we still want to transpose??? I think so...
    coords = properties.coords.T

    # NOTE: skip background corrections. Leave that for data analysis step.
    for ch, index in ch_to_index.items():
        # row['total_intensity_{}'.format(ch)] = stack[index, region_number, coords[0], coords[1]].sum()
        # TODO: make sure that the coords[0] & 1 didn't get Transposed! (when going from C to F ordering...)
        row["total_intensity_{}".format(ch)] = stack[
            coords[0], coords[1], region_number, index
        ].sum()

    row["info_filename"] = name
    row["info_fov"] = fov
    row["info_frame"] = frame
    row["info_z_level"] = z_level
    row["info_region_set_number"] = region_set_number  # e.g. a row of trenches
    row["info_region_number"] = region_number  # e.g. a trench
    row["info_label"] = properties.label
    row["info_seg_channel"] = sc

    row["geometry_Area"] = properties.area
    row["geometry_Orientation"] = properties.orientation
    row["geometry_Perimeter"] = properties.perimeter

    # FIXME very small numbers --> math domain error.
    # Shouldn't the skimage library handle this issue and just return NaN?
    # The values might drop below zero due to floating-point error,
    # and then a square root of a neg number causes the error.
    # Maybe we can wrap each of these in try clauses, and write NaN if they throw exceptions.
    row["axis_length_major"] = properties.major_axis_length
    row["axis_length_minor"] = properties.minor_axis_length

    row["centroid_row"] = properties.centroid[0]
    row["centroid_col"] = properties.centroid[1]

    row["bounding_box_min_row"] = properties.bbox[0]
    row["bounding_box_min_col"] = properties.bbox[1]
    row["bounding_box_max_row"] = properties.bbox[2]
    row["bounding_box_max_col"] = properties.bbox[3]

    # Append the properties information to the table
    row.append()


def make_cell_type(channels, seg_channels, file_names):
    """
    Define the column types for the PyTable: this stores segmented cell information
    NOTE: calling it a "Cell" because it stores information about biological cells,
    not because we're working with cells in a table.
    """
    # Define the PyTables column types using a dictionary
    column_types = {
        "info_fov": tables.UInt16Col(),
        "info_frame": tables.UInt16Col(),
        "info_z_level": tables.UInt16Col(),
        # e.g. the trench row number
        "info_region_set_number": tables.UInt16Col(),
        # e.g. trench number. The "region" within the whole image.
        "info_region_number": tables.UInt16Col(),
        # The labeled, connected component within the region
        "info_label": tables.UInt16Col(),
        "geometry_Area": tables.UInt16Col(),
        "geometry_Orientation": tables.Float32Col(),
        "geometry_Perimeter": tables.Float32Col(),
        # See note below about math domain errors.
        "axis_length_major": tables.Float32Col(),
        "axis_length_minor": tables.Float32Col(),
        "centroid_row": tables.UInt16Col(),
        "centroid_col": tables.UInt16Col(),
        "bounding_box_min_row": tables.UInt16Col(),
        "bounding_box_min_col": tables.UInt16Col(),
        "bounding_box_max_row": tables.UInt16Col(),
        "bounding_box_max_col": tables.UInt16Col(),
    }

    # NOTE if using flat field correction, make these floats
    # Otherwise, integers
    for c in channels:
        c = c.decode("utf-8")
        # Total intensities
        column_types["total_intensity_{}".format(c)] = tables.UInt32Col()

    # Need to know the maximum string length to make a string column
    max_channel_name_length = get_max_length(seg_channels)
    column_types["info_seg_channel"] = tables.StringCol(max_channel_name_length)

    max_filename_length = get_max_length(file_names)
    column_types["info_filename"] = tables.StringCol(max_filename_length)

    # Yeah, this syntax is really goofy
    return type("Cell", (tables.IsDescription,), column_types)


def get_max_length(items):
    """
    Input a list, return the length of the longest element
    """
    longest = 0

    for n in items:
        l = len(n)
        if l > longest:
            longest = l

    return longest


def subtract_background_from_coords(coords, image, background):
    """
    Input a list of pixel co-ordinates (as per skimage.measure.regionpprops),
    an associated image (2d numpy array) to read intensities from, and and a background signal value.
    Return the sum of the pixel intensities minus the product of the background value and the number of pixels.
    NOTE is it faster to use the coords[] array, or to "multiply" by the mask?
    """
    summed_intensity = image[coords[0], coords[1]].sum()
    total_background = len(coords) * background

    return summed_intensity - total_background


def subtract_background_from_region(image, image_mask, background):
    """
    Input an image, a binary mask, and a background value per pixel
    Return the sum of the pixel intensities minus the product of the background value and the number of pixels.
    """
    # Use mask to zero out unwanted pixels, and sum the intensities of the remaining ones
    summed_intensity = (image * image_mask).sum()

    # Calculate background signal
    # Sum of a mask is equivalent to the number of valid pixels, if the mask is as integer and not bool.
    total_background = image_mask.astype(numpy.uint8).sum() * background

    return summed_intensity - total_background


def merge_tables(in_file, out_file, channels, seg_channels, file_names):
    """
    Merge the tables at the very end
    """
    Cell = make_cell_type(channels, seg_channels, file_names)

    h5file_in = tables.open_file(in_file, mode="r")

    h5file_out = tables.open_file(out_file, mode="w")
    big_table = h5file_out.create_table("/", "concatenated_measurements", Cell)

    # NOTE: had to resort to much more complex code, and not a single walk_nodes,
    # because walk_nodes won't visit external links.
    ### See: https://stackoverflow.com/questions/24891014/append-all-rows-from-one-table-to-another-using-pytables
    # for small_table in h5file_in.walk_nodes(where="/", classname="Table"):
    for node1 in h5file_in.walk_nodes(where="/"):
        if isinstance(node1, tables.link.ExternalLink):
            node1 = node1()
            for node2 in h5file_in.walk_nodes(where=node1):
                if isinstance(node2, tables.link.ExternalLink):
                    node2 = node2()
                    for node3 in h5file_in.walk_nodes(where=node2):
                        if isinstance(node3, tables.link.ExternalLink):
                            node3 = node3()

                            for small_table in h5file_in.walk_nodes(
                                where=node3, classname="Table"
                            ):
                                for row in small_table.iterrows():
                                    big_table.append([row[:]])

    # Done copying trench properties
    big_table.flush()
    big_table.close()

    h5file_in.close()
    h5file_out.close()
