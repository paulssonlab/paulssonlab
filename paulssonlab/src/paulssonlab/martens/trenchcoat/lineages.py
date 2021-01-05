import numpy
import tables
import pandas
from dataframe_conversion import write_dataframe_to_hdf5

"""
Input cell properties for segmented cells and infer lineages based on several key assumptions:
- cells only exit through one end of a trench
- trenches are aligned vertically in the image. Cells move top to bottom,
  or bottom to top, not left to right or right to left.
- cells do not:
    - disappear
    - become smaller
    - change order
    - merge
- a division occurs when a cell becomes smaller than before by a minimum amount (defined by the user)
- the frequency between time frames is rapid enough, relative to the time
  between divisions, to not miss divisions
- Ignore z-levels (assume that all information is stored within a single z-level).

Also provide a few useful methods for re-labeling masks given the new lineages.

It's also important:
- not to have small (e.g. ~5 pixel) objects (e.g. disjointed cell fragments)
- not to have cells from neighboring trenches, as can happen if the images are
  allowed to extend beyond the trenches and cells can bend sideways.

Terminology:
- A kymograph is the list of lists of cell labels for a given trench.
  The first dimension is the time frame, and the second dimension is the labels within the frame.
- Each cell at time frame zero defines its own "lineage." This is because we can't know
  how these cells are interrelated prior to the beginning of the imaging.
- Each lineage has its own set of parent/progeny relationships.

Ultimate output is a new HDF5 table associating the relevant information
(FOV, frame, etc., and label, lineage_id, progeny_id).
Another program could then be used to merge this information with other data tables, as needed.
NOTE undecided if it isn't better instead to perform this merger now,
& to use pandas to write the merged HDF5 table.

TODO:
- check whether a cell is the last cell in a trench and maybe if its centroid or its
  tip furthest from the exit moves: this could help to avoid calling cell divisions
  when a cell is smaller because it's exiting the trench
- check whether two cells were accidentally merged due to a segmentation error.
  One way is to have a maximum doubling time, and if a doubling was called which is less
  than this doubling time, then re-merge the cells and re-calculate the kymograph & lineages.
  This could be particularly hairy if the merger happens between distinct lineages (and not
  just cells within a lineage).
"""


def make_kymograph_from_centroids(table, seg_channel, is_top_row):
    """
    1. Input dataframe slice corresponding to a single trench with a given File, FOV & trench number,
       a segmentation channel, and whether the row of trenches is top or bottom.
       Output a "kymograph" stored as a list of lists,
       and top_row (top_row == True) or bottom_row (top_row == False):
    a. Each sub-list represents the cells at each subsequent time frame.
    b. Its contents are the cells, ordered from the feeding channel to the end of the trench:
       re-order (label) them by comparing centroids (safe),
       or possibly assuming that their labeling was in this order (risky).
    c. A "cell" is stored in said list as the connected component region's label (integer).
    """
    kymograph = []

    # Get unique frame numbers from the table, in sorted order
    # NOTE this code assumes no gaps (indexing by lists), and in sorted frame order
    # TODO how to check for gaps? What to do if there are gaps?
    # One fix is to index with a dictionary.
    # Or, could we:
    # - pre-allocate an array given by the max frame number?
    # - use enumerate & explicit indexing instead of append
    # - or, still use append, but if there's a gap, append a None value,
    #   and then add code later to handle the Nones as needed.

    unique_frames = sorted(table["info_frame"].unique())

    # Iterate by time frame
    for f in unique_frames:
        # Grab the rows which belong to this frame by slicing the dataframe
        # Copy so that we can apply the sort
        cells_this_frame = table[table["info_frame"] == f].copy()

        # Sort the sliced dataframe by centroid row
        # NOTE: Ascending will be False (top) or True (bottom), depending on top or bottom row of trenches.
        cells_this_frame.sort_values(
            by=["centroid_row"], inplace=True, ascending=is_top_row
        )

        # Iterate in this new order and append the labels to a new list
        cell_order = []
        for c in cells_this_frame.itertuples():
            label = c.info_label
            length = c.axis_length_major
            this_cell = {"label": label, "cell_length": length}

            cell_order.append(this_cell)

        # And append this list to the kymograph
        kymograph.append(cell_order)

    return kymograph


def run_kymograph_to_lineages(
    kymograph, file, fov, seg_channel, lineage_row, trench_row, trench, length_buffer
):
    """
    Input labeled kymograph for a given FOV / trench, and the row for appending to the table.
    Write labeled cell lineage information to a table.
    This algorithm makes assumptions: see above.
    Allow length_buffer pixel margin of error +/- cell length b/w time frames,
    for determining whether or not a division took place.
    """
    # In the zeroth time frame, define cells as lineage progenitors.
    # For each, run the lineage algorithm on the remaining kymograph time frame entries.
    lineage = 0

    # Start at frame zero
    for cell in kymograph[0]:
        stack = []

        # Initialize the stack with this first cell (# 1).
        # Think of it as a "job" to be performed later.
        job = {
            "frame": 0,
            "lineage": lineage,
            "progeny_id": 1,
            "label": cell["label"],
            "cell_length": cell["cell_length"],
        }

        stack.append(job)

        # While it's not empty
        while stack:
            previous_cell = stack.pop()

            # First, write it to the table
            write_cell_to_lineage_table(
                previous_cell, lineage_row, file, fov, seg_channel, trench_row, trench
            )

            # Now look towards the next time frame
            this_progeny_id = previous_cell["progeny_id"]
            next_frame = previous_cell["frame"] + 1

            # NOTE do we check that there is a next frame here?
            # i.e. if (next_frame < len(kymograph):
            # If not, then there's no point in continuing, and no other cells to look for.

            # Are there more frames to process?
            # Are there any cells in the next frame?
            if (next_frame < len(kymograph)) and (kymograph[next_frame]):
                # Load the next frame from the kymograph
                next_frame_group = kymograph[next_frame]

                # Get the top cell: contains a label and a major axis length
                top_cell = next_frame_group.pop(0)

                # Is the top cell shorter than the cell in the previous frame? -> There was a cell division
                # NOTE if a cell is slowly sliding out of the trench, then
                # it will be mis-characterized as dividing, as it becomes shorter.
                # To fix this, can try to compare the position of the centroid or
                # the tip of the cell & see whether it moved?

                if (
                    top_cell["cell_length"] + length_buffer
                    < previous_cell["cell_length"]
                ):

                    # If the cell is ~ 2x larger than that in the previous frame, then
                    # there might be a cell merging issue (segmentation problem).
                    if top_cell["cell_length"] > 2 * previous_cell["cell_length"]:
                        print("Warning: possible cell merger detected.")

                    # Is there also a next (sister) cell, from this division? Add it to the stack first.
                    if next_frame_group:
                        # Contains a label and a major axis length
                        next_cell = next_frame_group.pop(0)

                        job = {
                            "frame": next_frame,
                            "lineage": lineage,
                            "progeny_id": this_progeny_id * 2 + 1,
                            "label": next_cell["label"],
                            "cell_length": next_cell["cell_length"],
                        }

                        stack.append(job)

                    # Now add the top cell to the stack
                    job = {
                        "frame": next_frame,
                        "lineage": lineage,
                        "progeny_id": this_progeny_id * 2,
                        "label": top_cell["label"],
                        "cell_length": top_cell["cell_length"],
                    }

                    stack.append(job)

                # Top cell is not shorter -> there was not a cell division
                # So don't change the progeny_id
                else:
                    job = {
                        "frame": next_frame,
                        "lineage": lineage,
                        "progeny_id": this_progeny_id,
                        "label": top_cell["label"],
                        "cell_length": top_cell["cell_length"],
                    }

                    stack.append(job)

            # No cells: do nothing!

        # End while loop: stack is empty, can move on to the next lineage, with a new stack and any remaining cells left in the kymograph.
        lineage += 1


def write_cell_to_lineage_table(
    cell, lineage_row, file, fov, seg_channel, trench_row, trench_number
):
    """
    Write a cell's lineage information to the table
    """
    lineage_row["info_file"] = file
    lineage_row["info_fov"] = fov
    lineage_row["info_seg_channel"] = seg_channel
    lineage_row["info_frame"] = cell["frame"]
    lineage_row["info_row_number"] = trench_row
    lineage_row["info_trench_number"] = trench_number
    lineage_row["lineage"] = cell["lineage"]
    lineage_row["progeny_id"] = cell["progeny_id"]
    lineage_row["info_label"] = cell["label"]

    lineage_row.append()


def make_lineage_table_type(filename_itemsize, segchannel_itemsize):
    """
    Helper function for writing tables to HDF5
    """
    # Define the PyTables column types using a dictionary
    column_types = {
        "info_file": tables.StringCol(filename_itemsize),
        "info_fov": tables.UInt16Col(),
        "info_frame": tables.UInt16Col(),
        "info_row_number": tables.UInt16Col(),
        "info_trench_number": tables.UInt16Col(),
        "info_seg_channel": tables.StringCol(segchannel_itemsize),
        "lineage": tables.UInt16Col(),
        "progeny_id": tables.UInt16Col(),
        "info_label": tables.UInt16Col(),
    }

    return type("Lineage", (tables.IsDescription,), column_types)


def find_ancestors(progeny_id):
    """
    Calculate the ancestor_ids for all ancestors of a given cell, and return them as a list.
    """
    ancestor_ids = []

    if progeny_id < 1:
        return None

    while progeny_id > 1:
        # Rely on integer division for odd-numbered ids to give even-numbered ancestors
        progeny_id //= 2
        ancestor_ids.append(progeny_id)

    return ancestor_ids


def find_descendents(
    desired_fov,
    desired_row,
    desired_trench,
    desired_lineage,
    desired_progeny_id,
    lineage_table,
):
    """
    Search the lineage table for all descendents of a given lineage & progeny_id
    TODO: make an alternative version which only returns the final members of the lineages.
    This could be useful for splitting the data up into sub-lineages in an intuitive way.
    """
    descendent_ids = []

    if progeny_id < 1:
        return None

    # Make a dataframe for easier searching, & just keep the rows for this particular fov, trench, lineage
    query = """(fov        == desired_fov)    & \
               (trench_row == desired_row)    & \
               (trench     == desired_trench) & \
               (lineage    == desired_lineage)"""

    matches = lineage_table.read_where(query)
    df = pandas.DataFrame(matches)

    # Use a stack to queue up searches for putative daughters
    stack = []

    putative_daughter_one = desired_progeny_id * 2
    putative_daughter_two = putative_daughter_one + 1

    stack.append(putative_daughter_one)
    stack.append(putative_daughter_two)

    # Exhaustively search all possible daughters
    while stack:
        next_descendent = stack.pop()

        # It exists!
        if next_descendent in df.progeny_id.values:
            # Add to confirmed list of descendents
            descendent_ids.append(next_descendent)

            # Add two putative daughters to the stack
            putative_daughter_one = next_descendent * 2
            putative_daughter_two = putative_daughter_one + 1

            stack.append(putative_daughter_one)
            stack.append(putative_daughter_two)

    return descendent_ids


def get_matching_cell_measurements_by_progeny_id(
    desired_lineage,
    progeny_ids,
    desired_fov,
    desired_trench,
    desired_row,
    lineage_table,
    measurements_table,
):
    """
    Input a list of lineage, progeny_ids & an FOV & trench, return their information from the cell_measurements table.
    Requires passing a query through the lineage_table first (using fov, trench, lineage, and then sub querying progeny_id),
    to associate id with the labels, then using those labels (and fov & trench) to grab out the specific cell table entries.
    This is helpful for e.g. grabbing all the cell info for a lineage, assuming that the ids represent cells from a same lineage.
    """
    # 1. Find all the progeny from this fov/trench/lineage
    matching_lineage = pandas.DataFrame(
        lineage_table.read_where(
            """(fov == desired_fov) & (trench_row == desired_row) & (trench == desired_trench) & (lineage == desired_lineage)"""
        )
    )
    # TODO: possible speedup for ints?
    # https://stackoverflow.com/questions/22485375/efficiently-select-rows-that-match-one-of-several-values-in-pandas-dataframe
    matching_progeny = matching_lineage[
        matching_lineage["progeny_id"].isin(progeny_ids)
    ]

    # 2. Find all the cell measurements which match fov & trench
    matching_measurements = pandas.DataFrame(
        measurements_table.read_where(
            """(info_fov == desired_fov) & (info_trench_row == desired_row) & (info_trench_number == desired_trench)"""
        )
    )

    # 3. For each matching cell from (1), associate its label info in (2) to get its measurements.
    # Presumably there should be exactly the same number of items in each query,
    # so the type of merge (inner, outer, left, right) won't make a difference?
    # Label alone won't work, because the labels are reset for every time frame...
    # but combine with 'frame' and it should work, right?
    merged = pandas.merge(
        matching_progeny,
        matching_measurements,
        left_on=["frame", "label"],
        right_on=["info_frame", "info_label"],
        how="inner",
    )

    # Return a Pandas dataframe out of the measurements_table
    return merged


def get_matching_labels(
    file, fov, frame, trench_row, trench, lineage, prog_id, dataframe
):
    """
    Given file, fov, frame, trench, lineage, progeny_id,
    and a dataframe which matches lineage info. with cell labels (from cell measurements),
    return the original (old) label from the mask, and the corresponding progeny_id.
    """
    matching_entry = dataframe[
        (dataframe["info_file"] == file)
        & (dataframe["info_fov"] == fov)
        & (dataframe["info_frame"] == frame)
        & (dataframe["info_trench_row"] == trench_row)
        & (dataframe["info_trench_number"] == trench)
        & (dataframe["lineage"] == lineage)
        & (dataframe["progeny_id"] == prog_id)
    ]

    if not matching_entry.empty:
        old_label = int(matching_entry.iloc[0].at["label"])
        new_label = int(matching_entry.iloc[0].at["progeny_id"])
        return (old_label, new_label)

    else:
        return None


def relabel_mask(mask, old_label, new_label, alternative_label):
    """
    Input a mask with labeled pixels.
    Set all pixels equal to the old label to the new label,
    and those not equal to the alternative label.
    This is a way to "highlight" a particular label.
    """
    new_mask = (mask == old_label) * new_label + (
        (mask != old_label) & (mask != 0)
    ) * alternative_label
    return new_mask


def relabel_mask_mono(mask, label):
    """
    Re-label all non-zero pixels with the given label
    """
    new_mask = (mask != 0) * label
    return new_mask


def relabel_mask_complete(mask, old_labels, new_labels):
    """
    For each label, replace it with a new label.
    Input a mask with labeled pixels.
    """
    new_mask = numpy.zeros(mask.shape, dtype=mask.dtype)
    for (o, n) in zip(old_labels, new_labels):
        new_mask += (mask == o) * n

    return new_mask


def get_matching_labels_global(
    file, fov, frame, trench_row, trench, lineage, dataframe
):
    """
    Given file, fov, frame, trench, lineage,
    and a dataframe which matches lineage info. with cell labels (from cell measurements),
    return the original (old) label from the mask, and ALL progeny_ids.
    """
    old_labels = []
    new_labels = []

    matching_entry = dataframe[
        (dataframe["info_file"] == file)
        & (dataframe["info_fov"] == fov)
        & (dataframe["info_frame"] == frame)
        & (dataframe["info_trench_row"] == trench_row)
        & (dataframe["info_trench_number"] == trench)
        & (dataframe["lineage"] == lineage)
    ]

    if not matching_entry.empty:
        for r in matching_entry.itertuples():
            old_label = r.label
            new_label = r.progeny_id

            old_labels.append(old_label)
            new_labels.append(new_label)

    else:
        return None

    # Not empty
    # FIXME is this check redundant?
    if old_labels:
        return (old_labels, new_labels)

    else:
        return None


### Main


def main_lineages_function(infile, outfile, length_buffer, seg_channel):
    """
    infile contains cell measurements table
    outfile contains a new table with the lineage information
    length_buffer for calling cell divisions, in pixels
    seg channel used to identify & label the cells TODO make this a list and process multiple seg channels?
    """
    # Open the measurements table as a Pandas dataframe
    # NOTE this code assumes that the data were processed
    # after re-numbering the trenches to account for stage drift over time.
    # (see: renumber_trenches.py)
    # df = pandas.read_hdf(infile, "/cells_in_labeled_trenches")
    h5_in = tables.open_file(infile, "r")
    # NOTE could also get the unique set of file names etc. from the metadata
    # in the original HDF5 file. This would avoid having to load the whole
    # data table into memory.
    table = h5_in.get_node("/cells_in_labeled_trenches")
    filename_itemsize = table.description.info_file.itemsize
    segchannel_itemsize = table.description.info_seg_channel.itemsize
    df = pandas.DataFrame(table.read())
    h5_in.close()

    # If there were multiple segmentation channels, just keep the one
    # TODO: here and elsewhere, use tables read_where() to more quickly weed out
    # the unwanted rows, potentially reducing the memory footprint.
    # Must compare against the bytes version of the string
    df = df[df["info_seg_channel"] == bytes(seg_channel, "utf-8")]

    # Define the table row type for writing lineage information
    Lineage = make_lineage_table_type(filename_itemsize, segchannel_itemsize)

    # Create the lineage information table
    h5_out = tables.open_file(outfile, "w")
    lineage_table = h5_out.create_table(
        "/",
        "lineages",
        Lineage,
        "Lineage_information",
        createparents=True,
        filters=tables.Filters(complevel=1, complib="zlib"),
    )

    # For appending entries
    lineage_row = lineage_table.row

    # TODO it might be possible to walk through the dataframe using
    # groupby and apply or map etc., but for now let's just loop
    # something like:
    # df.groupby(["info_file", "info_fov", "info_row_number", "corrected_trench_label"]).apply(lambda x: run_kymo_func(x))
    # def run_kymo_func(dataframe_slice): # call make_kymograph_from_centroids and run_kymograph_to_lineages
    # Note that the df isn't modified in the code, it just writes a new HDF5 file.
    # Alternatively, could add a new column to the original df, with the lineage_id and progeny_id,
    # however I have questions about the memory footprint vs. using PyTables.
    # Similarly, we are forced to load the whole table into memory above.
    # TODO Is there no way to selectively query portions of an on-disk HDF5 table written by Pandas?
    files = df["info_file"].unique()
    for file in files:
        this_file = df[(df["info_file"] == file)]
        fovs = this_file["info_fov"].unique()
        for fov in fovs:
            this_fov = this_file[(this_file["info_fov"] == fov)]
            trench_rows = this_fov["info_row_number"].unique()
            for row in trench_rows:
                this_row = this_fov[(this_fov["info_row_number"] == row)]
                trenches = this_row["corrected_trench_label"].unique()
                # If top row, set True, else False
                # This assumes an alternating even/odd numbering scheme.
                is_top_row = row % 2 == 0
                for trench in trenches:
                    this_trench = this_row[
                        (this_row["corrected_trench_label"] == trench)
                    ]

                    # Returns a list of lists for storing & popping cells in the same trench, but over time.
                    # Kymograph now contains cell labels & lengths over time
                    kymograph = make_kymograph_from_centroids(
                        this_trench, seg_channel, is_top_row
                    )

                    # Process the kymograph to calculate lineages & write them to a table
                    # NOTE "lineage_row" is the "row" in the data table where we write the results to HDF5.
                    # "row" is the number of the row of cells in the mother machine (0, 1, etc. for top, bottom...).
                    run_kymograph_to_lineages(
                        kymograph,
                        file,
                        fov,
                        seg_channel,
                        lineage_row,
                        row,
                        trench,
                        length_buffer,
                    )

    # Done writing lineage information
    # This table now contains all the information which can be used to map lineages back onto the kymograph.
    lineage_table.flush()
    h5_out.close()

    # Now, open the lineages table & create a new table merged with the cell measurements.
    # Write this to disk as well.
    h5_lineages = tables.open_file(outfile, "r")
    df_lineages = pandas.DataFrame(h5_lineages.get_node("/lineages").read())

    # Perform the merger
    merge_on = [
        "info_file",
        "info_fov",
        "info_row_number",
        "info_trench_number",
        "info_frame",
        "info_seg_channel",
        "info_label",
    ]
    merged = pandas.merge(df, df_lineages, on=merge_on)

    # Convert the columns to bytes ~ for whatever reason, this is required,
    # otherwise they are generic "objects" in the dataframe
    merged["info_file"] = merged["info_file"].astype(numpy.bytes_)
    merged["info_seg_channel"] = merged["info_seg_channel"].astype(numpy.bytes_)

    # TODO pass in another out file, OR pass in a directory and create custom outfiles?
    write_dataframe_to_hdf5(merged, "cells_and_lineages.h5", "cells_and_lineages")

    # Done
    h5_lineages.close()
