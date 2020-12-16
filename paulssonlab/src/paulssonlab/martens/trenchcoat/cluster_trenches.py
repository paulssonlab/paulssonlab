#!/usr/bin/env python

import tables
import pandas
import numpy
import sklearn.cluster as cluster
from metadata import get_metadata

"""
TODO:

- Put everything in a main function which can be invoked by Click
- Standardize on what serves as input ~ which files to open, params etc. Use a YAML file or not needed?
- Unit conversion from microns to pixels
- Move some code into subroutines?
- It seems to run very quickly, so don't worry about speed optimizations.

"""


def predict_clusters(x, column_name, eps, min_samples):
    """
    Only cluster using "min_row", which is the x-coordinate of the trenches.
    This will work particularly well if we can split rows of trenches apart
    thanks to the trench detection.
    Otherwise, it would be necessary to cluster in 2D spac using:
    numpy.column_stack((x["min_row"], x["min_col"]))

    Leaving this here in case it's useful
    # Use both the x & y info
    data = numpy.column_stack([x["min_row"], x["min_col"]])
    """
    # Reshape required because we use a single dimension for clustering
    data = numpy.array(x[column_name])
    data = data.reshape(-1, 1)

    ms = cluster.DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1)
    clusters = ms.fit_predict(data)

    x["corrected_trench_label"] = clusters

    return x


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


### Begin Main

# Load the FOV information
lane_number = 2
hdf5_file_path = "HDF5_lane_{}/data.h5".format(lane_number)
h5file = tables.open_file(hdf5_file_path, "r")
table_timestamps = h5file.get_node("/Metadata/File_run/fov_metadata")
ts_data = pandas.DataFrame(table_timestamps.read())

metadata_node = h5file.get_node("/Metadata/File_run")()
metadata = get_metadata(metadata_node)
h5file.close()

# Don't care about these columns
ts_data = ts_data.drop(["info_pfs_status", "info_pfs_offset"], axis=1)
normalize_stage_xyz_and_convert_to_pixels(ts_data, metadata["pixel_microns"])

# Load the trench regions information
infile_reg = "REG_lane_{}/regions_tables.h5".format(lane_number)
h5_reg = tables.open_file(infile_reg, "r")
reg_t = h5_reg.get_node("/trench_coords").read()
h5_reg.close()
df_reg_t = pandas.DataFrame(reg_t)

# Merge the trench regions coordinates table with the timestamps / xyz table
merged = pandas.merge(df_reg_t, ts_data, on=["info_fov", "info_frame"], how="left")

# Correct trench bounding boxes for stage drift
merged["min_row"] = merged["min_row"].astype(numpy.int32) - merged["info_x"]
merged["max_row"] = merged["max_row"].astype(numpy.int32) - merged["info_x"]
merged["min_col"] = merged["min_col"].astype(numpy.int32) + merged["info_y"]
merged["max_col"] = merged["max_col"].astype(numpy.int32) + merged["info_y"]

### Cluster the trenches using their min_row (max_row could work too)

# Distance between clusters, in pixels.
# If this value is too large, then clusters will be merged.
# If it is too small, then clusters will be split.
# The value therefore represents a compromise between residual
# drift between images, and the space between trenches.
eps = 10

# Minimum number of samples for a cluster to be considered not noise.
# This value depends on the size of the dataset insofar as that
# the number of frames determines the maximum number of members per cluster.
# Typically, a trench which is borderline will have fewer samples, and
# this value can be used to decide how often a trench needs to appear
# to makethe cut.
# Trenches which don't belong to large enough clusters are assigned group "-1".
min_samples = 100

# Result of the operation is stored as a new column in the original dataframe
# called "corrected_trench_label".
result = merged.groupby(
    ["info_file", "info_fov", "info_z_level", "info_row_number"]
).apply(lambda x: predict_clusters(x, "min_row", eps, min_samples))

# Load the cell measurements table
h5_seg = tables.open_file(
    "SEG_lane_{}/TABLES/tables_merged.h5".format(lane_number), "r"
)
properties_table = h5_seg.get_node("/concatenated_measurements")
df_cell_props = pandas.DataFrame(properties_table.read())

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
    result, df_cell_props, left_on=left_on, right_on=right_on
)

# These columns are now redundant
merge_tr_cell_total = merge_tr_cell_total.drop(
    ["info_filename", "info_region_number", "info_region_set_number"], axis=1
)

# Convert some of the column names. They now are more accessible, and can be written to HDF5 (ironically).
merge_tr_cell_total["info_file"] = merge_tr_cell_total["info_file"].map(
    lambda x: x.decode("utf-8")
)
merge_tr_cell_total["info_seg_channel"] = merge_tr_cell_total["info_seg_channel"].map(
    lambda x: x.decode("utf-8")
)

# Write to disk
outfile = "test_hdf5_datatable.h5"
merge_tr_cell_total.to_hdf(
    outfile, key="cells_in_labeled_trenches", mode="w", complevel=1, format="table"
)
