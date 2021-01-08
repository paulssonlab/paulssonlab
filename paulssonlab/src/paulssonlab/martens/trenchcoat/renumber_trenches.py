#!/usr/bin/env python

import tables
import pandas
import numpy
from metadata import get_metadata
from dataframe_conversion import write_dataframe_to_hdf5

"""
TODO:
- Unit conversion from microns to pixels?
"""


def adjust_trench_number(df, threshold):
    """
    FIXME This still isn't quite right (fails sometimes). Is it a matter of tweaking the threshold, or a deeper problem?

    threshold: how much horizontal shift (in pixels) before we declare that we gained or lost a trench.
    This parameter is key for determining if a shift occurred, and the validity of this
    algorithm ultimately depends on the accuracy of correcting for stage drift
    (e.g. using the drift recorded by the stage), and how close together the trenches are.
    Assumes that the shift is a most 1 trench width's worth. Unclear if or how well
    this would extend to larger drifts.
    """
    # For each file / fov / frame / z / row, obtain the absolute minimum min_row value.
    smallest_min_rows = df.groupby(
        ["info_file", "info_fov", "info_z_level", "info_row_number"]
    )["min_row"].min()
    # Convert from a series to a dataframe
    smallest_min_rows = smallest_min_rows.reset_index()

    # Go through the minimum min_row values
    for i in smallest_min_rows.itertuples():
        # Pick out the frames in the same FOV etc.
        this_df = df[
            (df["info_file"] == i.info_file)
            & (df["info_fov"] == i.info_fov)
            & (df["info_z_level"] == i.info_z_level)
            & (df["info_row_number"] == i.info_row_number)
        ]

        # In principle, using unique() allows for gaps in the frame numbering.
        all_frames = this_df["info_frame"].unique()
        for fr in all_frames:
            this_fr = this_df[this_df["info_frame"] == fr]
            # This is the very first min_row, belonging to trench 0 according to this FOV.
            zeroth_min_row = this_fr[this_fr["info_trench_number"] == 0][
                "min_row"
            ].iloc[0]

            # Is this the true trench 0, or has there been a shift?
            # Let's compare this min_row to the absolute smallest min_row,
            # and see if the difference falls within a threshold.
            diff = zeroth_min_row - i.min_row

            # If it does, then the difference is too small to signify a shift.
            # FIXME why are the values in the dataframe stored as floating point??
            # -> Because when it adds the first values, it must allocate an empty column,
            # and the remaining values must be assigned, so it assigns them NaN and forces
            # conversion to f64. Options are:
            # 1. converting back to int afterwards, or
            # 2. using a groupby method, which might preserve data types.
            if diff < threshold:
                df.loc[this_fr.index, ("corrected_trench_label")] = this_fr[
                    "info_trench_number"
                ]
            # Otherwise, we estimate that there is too large a gap, and so we increment all min_row
            # values by one throughout this time frame.
            else:
                df.loc[this_fr.index, ("corrected_trench_label")] = (
                    this_fr["info_trench_number"] + 1
                )

    # Returns the original dataframe with a new column, "corrected_trench_label".
    # Cast to uint16 TODO not have to cast
    df = df.astype({"corrected_trench_label": "uint16"})

    # And cast the file name to bytes to support writing to hdf5 without pandas
    # NOTE WTF why doesn't this work using the dictionary style above?
    df["info_file"] = df["info_file"].astype(numpy.bytes_)

    return df


def normalize_stage_xyz(ts_data):
    """
    Input a Pandas dataframe containing columns for fov, frame, x, y, z etc.
    *Modify* the dataframe with the x, y, z normalized by zeroth time frame
    on a per-FOV basis.

    - Group by FOV number
    - Subtract zeroth x,y,z element (first time frame) from all x,y,z elements
    - Convert series back to a dataframe and drop two indexing columns
    - Grab the x,y,z columns
    - Assign the result into the x,y,z positions in the original dataframe
    """
    ts_data.loc[:, ("info_x", "info_y", "info_z")] = (
        ts_data.groupby("info_fov")
        .apply(
            lambda a: a.loc[:, ("info_x", "info_y", "info_z")]
            - a.loc[:, ("info_x", "info_y", "info_z")].iloc[0]
        )
        .reset_index(0, 1)
        .loc[:, ("info_x", "info_y", "info_z")]
    )


def subtract_and_divide(a, px_mu):
    """
    Helper function for normalizing stage xyz values (see below)
    """
    # Subtract the zeroth x, y, z information from all rows
    diff = a.loc[:, ("info_x", "info_y", "info_z")].astype(numpy.int32) - a.loc[
        :, ("info_x", "info_y", "info_z")
    ].iloc[0].astype(numpy.int32)

    # Convert from microns to pixels
    diff /= px_mu
    diff = diff.astype(numpy.int16)

    return diff


def normalize_stage_xyz_and_convert_to_pixels(ts_data, px_mu):
    """
    Input a Pandas dataframe containing columns for fov, frame, x, y, z etc.
    *Modify* the dataframe with the x, y, z normalized by zeroth time frame
    on a per-FOV basis, and then converted to pixels from microns.

    - Group by FOV number
    - Subtract zeroth x,y,z element (first time frame) from all x,y,z elements
    - Convert series back to a dataframe and drop two indexing columns
    - Grab the x,y,z columns
    - Assign the result into the x,y,z positions in the original dataframe
    """
    ts_data.loc[:, ("info_x", "info_y", "info_z")] = (
        ts_data.groupby("info_fov")
        .apply(lambda a: subtract_and_divide(a, px_mu))
        .reset_index(0, 1)
        .loc[:, ("info_x", "info_y", "info_z")]
    )


def update_cell_measurements(
    adjusted_regions, cell_measurements_file, cell_measurements_outfile
):
    """
    Add the new trench coordinates to the cell measurements table.
    Also add the timestamps (they're handy).
    """
    # Load the cell measurements table
    h5_seg = tables.open_file(cell_measurements_file, "r")
    properties_table = h5_seg.get_node("/concatenated_measurements")
    df_cell_props = pandas.DataFrame(properties_table.read())
    h5_seg.close()

    # Merge the cell measurements with the re-labeled trenches
    left_on = [
        "info_file",
        "info_fov",
        "info_z_level",
        "info_row_number",
        "info_trench_number",
        "info_frame",
    ]
    right_on = [
        "info_filename",
        "info_fov",
        "info_z_level",
        "info_region_set_number",
        "info_region_number",
        "info_frame",
    ]
    merge_tr_cell_total = pandas.merge(
        adjusted_regions, df_cell_props, left_on=left_on, right_on=right_on
    )

    # These columns are now redundant
    merge_tr_cell_total = merge_tr_cell_total.drop(
        ["info_filename", "info_region_number", "info_region_set_number"], axis=1
    )

    # And we don't care about these: the regions file will always be more complete, as we did an inner merge.
    merge_tr_cell_total = merge_tr_cell_total.drop(
        ["min_row", "min_col", "max_row", "max_col", "info_x", "info_y", "info_z"],
        axis=1,
    )

    # Convert the info_file column back to numpy.bytes_ type
    merge_tr_cell_total["info_file"] = merge_tr_cell_total["info_file"].astype(
        numpy.bytes_
    )

    # And the same for the seg channel
    merge_tr_cell_total["info_seg_channel"] = merge_tr_cell_total[
        "info_seg_channel"
    ].astype(numpy.bytes_)

    # Write to disk
    write_dataframe_to_hdf5(
        merge_tr_cell_total, cell_measurements_outfile, "cells_in_labeled_trenches"
    )


def load_xyz_data(hdf5_file_path):
    """
    Load xyz & timestamps, normalize to zeroth timepoint.
    """
    h5file = tables.open_file(hdf5_file_path, "r")
    table_timestamps = h5file.get_node("/Metadata/File_run/fov_metadata")
    xyz_data = pandas.DataFrame(table_timestamps.read())
    metadata_node = h5file.get_node("/Metadata/File_run")()
    metadata = get_metadata(metadata_node)
    h5file.close()

    # Don't care about these columns
    xyz_data = xyz_data.drop(["info_pfs_status", "info_pfs_offset"], axis=1)

    # Normalize
    normalize_stage_xyz_and_convert_to_pixels(xyz_data, metadata["pixel_microns"])

    return xyz_data


def load_regions(regions_file):
    """
    Load trench regions into a dataframe.
    Fix the utf-8 weirdness in the info_file column.
    """
    # Load the trench regions information
    h5_reg = tables.open_file(regions_file, "r")
    reg = h5_reg.get_node("/trench_coords").read()
    h5_reg.close()
    df_regions = pandas.DataFrame(reg)

    # Convert some of the column names. They now are more accessible, and can be written to HDF5 (ironically).
    df_regions["info_file"] = df_regions["info_file"].map(lambda x: x.decode("utf-8"))

    return df_regions


def main_renumbering_function(
    in_file,
    regions_file,
    regions_outfile,
    cell_measurements_file,
    cell_measurements_outfile,
    threshold,
):
    """
    Use XYZ information in conjunction with trench coordinates to determine when
    a shift in trench numbering occurred.
    Also add the new labels to the cell measurements table.
    This function does not take into account the possibility that the trench detection
    procedure itself may be error-prone.
    """
    # Load the FOV information
    xyz_data = load_xyz_data(in_file)

    # Load the regions information
    df_regions = load_regions(regions_file)

    # Merge the trench regions coordinates table with the timestamps / xyz table
    regions = pandas.merge(
        df_regions, xyz_data, on=["info_fov", "info_frame"], how="left"
    )

    # ### KEY FUNCTION ###
    # Correct trench bounding boxes for stage drift
    # NOTE if we pass in actual images, then it may be possible to supplelement
    # this step with e.g. phase cross correlation or template matching.
    regions["min_row"] = regions["min_row"].astype(numpy.int32) - regions[
        "info_x"
    ].astype(numpy.int32)
    regions["max_row"] = regions["max_row"].astype(numpy.int32) - regions[
        "info_x"
    ].astype(numpy.int32)
    regions["min_col"] = regions["min_col"].astype(numpy.int32) + regions[
        "info_y"
    ].astype(numpy.int32)
    regions["max_col"] = regions["max_col"].astype(numpy.int32) + regions[
        "info_y"
    ].astype(numpy.int32)

    # ### KEY FUNCTION ###
    # Given the stage drift-corrected regions, determine for every FOV & Frame
    # whether or not a trench appeared or disappeared at the leftmost edge,
    # causing an offset by 1 in the trench numbering.
    adjusted_regions = adjust_trench_number(regions, threshold)

    # Write the newly labeled regions to disk
    write_dataframe_to_hdf5(adjusted_regions, regions_outfile, "regions_table")

    ### Add the new information to the cell measurements table
    if cell_measurements_outfile:
        update_cell_measurements(
            adjusted_regions, cell_measurements_file, cell_measurements_outfile
        )
