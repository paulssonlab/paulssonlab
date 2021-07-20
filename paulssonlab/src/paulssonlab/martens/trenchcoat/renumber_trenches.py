#!/usr/bin/env python
import os
import tables
import pandas
import numpy
from metadata import get_metadata, get_files_list, get_attribute_list
from dataframe_conversion import write_dataframe_to_hdf5

from multiprocessing import Pool
from skimage.registration import phase_cross_correlation

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
            # 3. store the result in a numpy array, and then add the numpy array as a new column at the end?
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
    diff = (
        a.loc[:, ("info_x", "info_y", "info_z")]
        - a.loc[:, ("info_x", "info_y", "info_z")].iloc[0]
    )

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
        ts_data.groupby(["info_file", "info_fov"])
        .apply(lambda a: subtract_and_divide(a, px_mu))
        .reset_index(0, 1)
        .loc[:, ("info_x", "info_y", "info_z")]
    )

    return ts_data


def update_cell_measurements(
    adjusted_regions, cell_measurements_file, cell_measurements_outfile
):
    """
    Add the new trench coordinates to the cell measurements table.
    Also add the timestamps (they're handy).
    """
    # Load the cell measurements table
    h5_seg = tables.open_file(
        os.path.join(cell_measurements_file, "TABLES/tables_merged.h5"), "r"
    )
    properties_table = h5_seg.get_node("/cell_measurements")
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
        "info_file",
        "info_fov",
        "info_z_level",
        "info_row_number",
        "info_trench_number",
        "info_frame",
    ]
    merge_tr_cell_total = pandas.merge(
        adjusted_regions, df_cell_props, left_on=left_on, right_on=right_on
    )

    # These columns are now redundant
    # TODO check whether the duplicate column is dropped automatically
    # merge_tr_cell_total = merge_tr_cell_total.drop(["info_file", "info_trench_number", "info_row_number"], axis=1)

    # And we don't care about these: the regions file will always be more complete, as we did an inner merge.
    merge_tr_cell_total = merge_tr_cell_total.drop(
        ["min_row", "min_col", "max_row", "max_col", "info_x", "info_y"], axis=1
    )

    # And depending on the method, maybe also drop info_z
    # FIXME this doesn't work?
    try:
        merge_tr_cell_total = merge_tr_cell_total.drop(["info_z"], axis=1)
    except:
        pass

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
        merge_tr_cell_total, cell_measurements_outfile, "cell_measurements"
    )

    # And then add the stored seg. params information to the new file, too
    print("Copying YAML parameters...")
    h5file = tables.open_file(
        os.path.join(cell_measurements_file, "MASKS/masks.h5"), "r"
    )
    seg_params_node = h5file.get_node("/Parameters", "seg_params.yaml")
    seg_params_str = seg_params_node.read()
    h5file.close()

    h5file_masks = tables.open_file(cell_measurements_outfile, "r+")
    h5file_masks.create_array(
        "/Parameters",
        "seg_params.yaml",
        obj=numpy.array(seg_params_str),
        createparents=True,
    )
    h5file_masks.close()


def load_xyz_data(hdf5_file_path):
    """
    Load xyz & timestamps, normalize to zeroth timepoint.
    """
    all_dfs = []
    h5file = tables.open_file(hdf5_file_path, "r")
    metadata_node = h5file.get_node("/Metadata")
    for n in h5file.list_nodes(metadata_node()):
        # Load the metadata for pixel microns
        metadata = get_metadata(n())
        table_timestamps = h5file.get_node(n, "fov_metadata")

        # Read into a dataframe
        xyz_data = pandas.DataFrame(table_timestamps.read())

        # Add a new column, which is the file name
        xyz_data["info_file"] = n._v_name

        # Don't care about these columns
        xyz_data = xyz_data.drop(["info_pfs_status", "info_pfs_offset"], axis=1)

        # Normalize
        xyz_data = normalize_stage_xyz_and_convert_to_pixels(
            xyz_data, metadata["pixel_microns"]
        )

        # Append to the list
        all_dfs.append(xyz_data)

    h5file.close()

    # Merge all the dataframes, from each file
    all_dfs = pandas.concat(all_dfs)

    return all_dfs


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


def cross_corr_offsets(file, fov, z_level, channel, df, in_file):
    """
    Input a dataframe with region information,
    a channel for comparing between frames,
    and the path to the H5 file which contains the images.

    Return X/Y offsets for every frame, relative to the zeroth
    (first available) frame.

    TODO is there a way to supplement the calculation with stage offsets,
    to give a little bit of help? (In case there is a very large offset
    which was mostly recorded by the stage).
    Can't figure out the proper way to crop the images accordingly before comparing.
    """
    h5_in = tables.open_file(in_file, "r")

    # Store all shifts here. At the end, will convert to Pandas dataframe.
    total_size = len(df)
    result = {
        "info_file": [],  # Storing bytes is trickier... use a regular Python list.
        "info_fov": numpy.empty(shape=total_size, dtype=numpy.uint16),
        "info_frame": numpy.empty(shape=total_size, dtype=numpy.uint16),
        "info_z_level": numpy.empty(shape=total_size, dtype=numpy.uint16),
        "info_x": numpy.empty(shape=total_size, dtype=numpy.int16),
        "info_y": numpy.empty(shape=total_size, dtype=numpy.int16),
    }

    # Get the unique values and sort.
    # Try to avoid issues if there are gaps in the frame numbers.
    # FIXME just sort the dataframe and walk through using itertuples
    all_frames = df["info_frame"].unique()
    all_frames = sorted(all_frames)

    # Pre-initialize, before the loop, with the first two frames.
    # This is because otherwise we'll run out of bounds on the
    # final frame, trying to look ahead into the eternal void.
    result["info_file"].append(file)
    result["info_fov"][0] = fov
    result["info_frame"][0] = all_frames[0]
    result["info_z_level"][0] = z_level
    result["info_x"][0] = 0
    result["info_y"][0] = 0

    # The initial frame.
    path = "/Images/{}/FOV_{}/Frame_{}/Z_{}/{}".format(
        file, fov, all_frames[0], z_level, channel
    )
    frame_b = h5_in.get_node(path).read()

    # Zeroth frame is special. Remove it.
    all_frames.pop(0)

    # Iterate the actual frame values, rather than their positions,
    # in case there are gaps between frames.
    # Always use the closest available frame for comparing.
    for i, fr in enumerate(all_frames):
        # Now, frame_b is reassigned as frame_a
        frame_a = frame_b

        # And we load a new frame_b
        path = "/Images/{}/FOV_{}/Frame_{}/Z_{}/{}".format(
            file, fov, fr, z_level, channel
        )
        frame_b = h5_in.get_node(path).read()

        # Compute the difference between the two images
        # NOTE frame_a comes first, or second?
        (shift_y, shift_x), error, diffphase = phase_cross_correlation(frame_a, frame_b)

        # Normalize this to the previous frame.
        # This ultimately has the effect of normalizing all shifts
        # to the zeroth frame.
        # Note that i enumerates from zero, but result already has
        # 1 element. So compare against element i, and add to element
        shift_x += result["info_x"][i]
        shift_y += result["info_y"][i]

        # Store the result
        result["info_file"].append(file)
        result["info_fov"][i + 1] = fov
        result["info_frame"][i + 1] = fr
        result["info_z_level"][i + 1] = z_level
        result["info_x"][i + 1] = shift_x
        result["info_y"][i + 1] = shift_y

    # Done!
    h5_in.close()
    return pandas.DataFrame.from_dict(result)


def start_cross_corr(in_file, channel):
    """
    Use the cross-correlation computed between pairs of images to determine the amount of shift.
    """
    files = get_files_list(in_file)

    h5_in = tables.open_file(in_file, "r")
    file_to_z = get_attribute_list(h5_in, "z_levels")
    h5_in.close()

    # Store individual DF results in here,
    # concatenate into a single DF at the end.
    all_offsets = []

    def log_result(result):
        # This is called whenever Pool returns a result.
        # Result_list is modified only by the main process, not the pool workers.
        # https://izziswift.com/multiprocessing-pool-when-to-use-apply-apply_async-or-map/
        all_offsets.append(result)

    # TODO figure out multi-processing
    with Pool() as pool:
        for file in files:
            h5_in = tables.open_file(in_file, "r")
            path = "/Metadata/{}/fov_metadata".format(file)
            table_fov_metadata = h5_in.get_node(path).read()
            h5_in.close()
            df_fov_metadata = pandas.DataFrame(table_fov_metadata)
            fovs = df_fov_metadata["info_fov"].unique()

            for fov in fovs:
                print("FOV {}".format(fov))
                df_this_fov = df_fov_metadata[df_fov_metadata["info_fov"] == fov]
                z_levels = file_to_z[file]

                for z_level in z_levels:
                    # Z-level information isn't stored in the regions
                    # table, but that doesn't have to stop us from
                    # calculating drift for every z-level.
                    args = [file, fov, z_level, channel, df_this_fov, in_file]
                    # offset = cross_corr_offsets(*args)
                    # all_offsets.append(offset)
                    pool.apply_async(cross_corr_offsets, args, callback=log_result)

        pool.close()
        pool.join()

        # Concatenate all the individual dataframe results into a single dataframe
        big_df = pandas.concat(all_offsets)

    return big_df


def main_renumbering_function(
    in_file,
    regions_file,
    regions_outfile,
    cell_measurements_file,
    cell_measurements_outfile,
    threshold,
    method,
    channel,
):
    """
    Use XYZ information in conjunction with trench coordinates to determine when
    a shift in trench numbering occurred.
    Also add the new labels to the cell measurements table.
    This function does not take into account the possibility that the trench detection
    procedure itself may be error-prone.
    """
    # Load the regions information
    df_regions = load_regions(regions_file)

    in_file = os.path.join(in_file, "data.h5")

    # Method 1
    # Use images to calculate offsets
    if method == "image":
        print("Using image cross correlations to calculate drift.")
        df_offsets = start_cross_corr(in_file, channel)
        # DEBUG
        df_offsets.to_hdf("image_result.h5", "data")

        # Merge the trench regions coordinates table with the computed img offsets
        regions = pandas.merge(
            df_regions,
            df_offsets,
            on=["info_file", "info_fov", "info_frame", "info_z_level"],
            how="left",
        )

    # Method 2
    # Load the FOV information
    elif method == "stage":
        print("Using stage information to calculate drift.")
        xyz_data = load_xyz_data(in_file)
        # DEBUG
        # xyz_data.to_hdf("stage_result.h5", "data")

        # Merge the trench regions coordinates table with the timestamps / xyz table
        regions = pandas.merge(
            df_regions, xyz_data, on=["info_file", "info_fov", "info_frame"], how="left"
        )
    else:
        print("Error: no valid method provided. Must be image or stage.")
        return

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
