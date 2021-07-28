#!/usr/bin/env python3

import numpy
import tables

"""
Functions for handling ND2 cell properties (measurements), including HDF5 tables.
"""


def init_properties_dict(
    channels, has_regions, size, max_len_filenames, max_len_seg_channels
):
    """
    Create a dict of empty lists for storing single cell measurements (properties).
    Eventually, this will be written to an HDF5 table.

    NOTE Would dynamic numpy arrays be helpful?
    https://github.com/maciejkula/dynarray

    Or, pass in size for a given batch:
    all the regions & masks for this particular task.

    Then, we collect a list of dicts of numpy arrays.
    Finally, we can either concatenate them into a dataframe, & write to disk,
    or just iterate them all.

    Note that it will be necessary to collect multiple such dicts of there are multiple rows,
    because these will be passed back to the same callback function.

    Considerations include the time it takes to convert to a dataframe,
    and the faster (vectorized) iterations that pandas can do compared
    to python for loops.

    The hope is that, by specifying the dtypes up front, that we will
    minimize conversion costs.

    FIXME how to determine the dtype for the 2 string type columns:
    1. info_seg_channel
    2. info_file

    In principle, we know the max lengths in advance,
    and could pass these all the way down through the functions.
    Should it be numpy.bytes_, or ???
    """
    results = {
        "info_file": numpy.empty(size, dtype="S{}".format(max_len_filenames)),
        "info_fov": numpy.empty(size, dtype=numpy.uint16),
        "info_frame": numpy.empty(size, dtype=numpy.uint16),
        "info_z_level": numpy.empty(size, dtype=numpy.uint16),
        "info_row_number": numpy.empty(size, dtype=numpy.uint16),
        "info_trench_number": numpy.empty(size, dtype=numpy.uint16),
        "info_label": numpy.empty(size, dtype=numpy.uint16),
        "info_seg_channel": numpy.empty(size, dtype="S{}".format(max_len_seg_channels)),
        "geometry_Area": numpy.empty(size, dtype=numpy.uint16),
        "geometry_Orientation": numpy.empty(size, dtype=numpy.float32),
        "geometry_Perimeter": numpy.empty(size, dtype=numpy.float32),
        "axis_length_major": numpy.empty(size, dtype=numpy.float32),
        "axis_length_minor": numpy.empty(size, dtype=numpy.float32),
        "centroid_row": numpy.empty(size, dtype=numpy.uint16),
        "centroid_col": numpy.empty(size, dtype=numpy.uint16),
        "bounding_box_min_row": numpy.empty(size, dtype=numpy.uint16),
        "bounding_box_min_col": numpy.empty(size, dtype=numpy.uint16),
        "bounding_box_max_row": numpy.empty(size, dtype=numpy.uint16),
        "bounding_box_max_col": numpy.empty(size, dtype=numpy.uint16),
    }

    for c in channels:
        results["total_intensity_{}".format(c)] = numpy.empty(size, dtype=numpy.float64)

    if has_regions:
        results["geometry_Max_Width_Area"] = numpy.empty(size, dtype=numpy.uint32)
        for c in channels:
            results["width_intensity_{}".format(c)] = numpy.empty(
                size, dtype=numpy.float64
            )

    return results


def add_properties(
    results,
    mask,
    file,
    fov,
    frame,
    z_level,
    row_number,
    trench_number,
    seg_channel,
    stack,
    ch_to_index,
    cell_properties,
    index,
):
    """
    Add a single cell's properties (from skimage measure) to the dict of lists.
    Modify the dictionary in-place (don't need to return or re-assign).
    """
    # NOTE skimage is changing their coordinates -- do we still want to transpose??? I think so...
    coords = cell_properties.coords.T

    # NOTE: skip background corrections. Leave that for data analysis step.
    for ch, i in ch_to_index.items():
        results["total_intensity_{}".format(ch)][index] = stack[
            coords[0], coords[1], trench_number, i
        ].sum()

    # Are there multiple regions?
    if stack.shape[2] > 1:
        # Filter out pixels that either belong to this cell
        # or to no cell at all (exclude other cells).
        not_other_cell = (
            mask[..., cell_properties.bbox[1] : cell_properties.bbox[3]]
            == cell_properties.label
        ) + (mask[..., cell_properties.bbox[1] : cell_properties.bbox[3]] == 0)

        # FIXME but this doesn't take into account the no_other_cell bit
        # Instead, make it the sum of no_other_cell?
        # row["geometry_Max_Width_Area"] = (stack.shape[0]) * (cell_properties.bbox[3] - cell_properties.bbox[1])
        results["geometry_Max_Width_Area"][index] = not_other_cell.sum()

        for ch, i in ch_to_index.items():
            rect_region = stack[
                ..., cell_properties.bbox[1] : cell_properties.bbox[3], trench_number, i
            ]
            rect_region *= not_other_cell
            results["width_intensity_{}".format(ch)][index] = rect_region.sum()

    results["info_file"][index] = file
    results["info_fov"][index] = fov
    results["info_frame"][index] = frame
    results["info_z_level"][index] = z_level
    results["info_row_number"][index] = row_number
    results["info_trench_number"][index] = trench_number
    results["info_label"][index] = cell_properties.label
    results["info_seg_channel"][index] = seg_channel

    results["geometry_Area"][index] = cell_properties.area
    results["geometry_Orientation"][index] = cell_properties.orientation
    results["geometry_Perimeter"][index] = cell_properties.perimeter

    # FIXME very small numbers --> math domain error.
    # Shouldn't the skimage library handle this issue and just return NaN?
    # The values might drop below zero due to floating-point error,
    # and then a square root of a neg number causes the error.
    # Maybe we can wrap each of these in try clauses, and write NaN if they thresults exceptions.
    results["axis_length_major"][index] = cell_properties.major_axis_length
    results["axis_length_minor"][index] = cell_properties.minor_axis_length

    results["centroid_row"][index] = cell_properties.centroid[0]
    results["centroid_col"][index] = cell_properties.centroid[1]

    results["bounding_box_min_row"][index] = cell_properties.bbox[0]
    results["bounding_box_min_col"][index] = cell_properties.bbox[1]
    results["bounding_box_max_row"][index] = cell_properties.bbox[2]
    results["bounding_box_max_col"][index] = cell_properties.bbox[3]


def write_properties_to_table_from_df(cell, table_row, columns):
    """
    Input a row from a dataframe, and an HDF5 table row element.
    Write all the values from one to the other.
    """
    for column in columns:
        table_row[column] = getattr(cell, column)

    table_row.append()


def write_properties_to_table(
    name,
    fov,
    frame,
    z_level,
    row_number,
    trench_number,
    properties,
    sc,
    row,
    stack,
    ch_to_index,
    mask,
):
    """
    Once a given trench has a segmentation mask, then write the properties of the masked region to the table.
    NOTE future versions of skimage will allow spitting out a properties object all at once, rather than lazily calculating them one at a time.
    """
    # NOTE skimage is changing their coordinates -- do we still want to transpose??? I think so...
    coords = properties.coords.T

    # NOTE: skip background corrections. Leave that for data analysis step.
    for ch, index in ch_to_index.items():
        # row['total_intensity_{}'.format(ch)] = stack[index, trench_number, coords[0], coords[1]].sum()
        # TODO: make sure that the coords[0] & 1 didn't get Transposed! (when going from C to F ordering...)
        row["total_intensity_{}".format(ch)] = stack[
            coords[0], coords[1], trench_number, index
        ].sum()

    # Are there multiple regions?
    if stack.shape[2] > 1:
        # Filter out pixels that either belong to this cell
        # or to no cell at all (exclude other cells).
        not_other_cell = (
            mask[..., properties.bbox[1] : properties.bbox[3]] == properties.label
        ) + (mask[..., properties.bbox[1] : properties.bbox[3]] == 0)

        # FIXME but this doesn't take into account the no_other_cell bit
        # Instead, make it the sum of no_other_cell?
        # row["geometry_Max_Width_Area"] = (stack.shape[0]) * (properties.bbox[3] - properties.bbox[1])
        row["geometry_Max_Width_Area"] = not_other_cell.sum()

        for ch, index in ch_to_index.items():
            rect_region = stack[
                ..., properties.bbox[1] : properties.bbox[3], trench_number, index
            ]
            rect_region *= not_other_cell
            row["width_intensity_{}".format(ch)] = rect_region.sum()

    row["info_file"] = name
    row["info_fov"] = fov
    row["info_frame"] = frame
    row["info_z_level"] = z_level
    row["info_row_number"] = row_number
    row["info_trench_number"] = trench_number
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


def make_cell_type(channels, seg_channels, file_names, has_regions):
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
        "info_row_number": tables.UInt16Col(),
        # e.g. trench number. The "region" within the whole image.
        "info_trench_number": tables.UInt16Col(),
        # The labeled, connected component within the region
        "info_label": tables.UInt16Col(),
        "geometry_Area": tables.UInt32Col(),
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

    for c in channels:
        # Total intensities
        column_types["total_intensity_{}".format(c)] = tables.Float64Col()

    # Add additional columns for mother machine data
    if has_regions:
        column_types["geometry_Max_Width_Area"] = tables.UInt32Col()

        for c in channels:
            # Total intensities spanning the entire width of the trench
            column_types["width_intensity_{}".format(c)] = tables.Float64Col()

    # Need to know the maximum string length to make a string column
    max_channel_name_length = get_max_length(seg_channels)
    column_types["info_seg_channel"] = tables.StringCol(max_channel_name_length)

    max_filename_length = get_max_length(file_names)
    column_types["info_file"] = tables.StringCol(max_filename_length)

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


def merge_tables(in_file, out_file, channels, seg_channels, file_names, has_regions):
    """
    Merge the tables at the very end
    FIXME something wrong with the final node3() de-referencing, having to do with File_{}/ component to the path?
    """
    Cell = make_cell_type(channels, seg_channels, file_names, has_regions)

    h5file_in = tables.open_file(in_file, mode="r")

    h5file_out = tables.open_file(out_file, mode="w")
    big_table = h5file_out.create_table(
        "/",
        "cell_measurements",
        Cell,
        filters=tables.Filters(complevel=1, complib="zlib"),
    )

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
