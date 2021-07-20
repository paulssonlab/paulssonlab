#!/usr/bin/env python
import tables
import pandas
import numpy
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
    lineage_row["corrected_trench_label"] = trench_number
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
        "corrected_trench_label": tables.UInt16Col(),
        "info_seg_channel": tables.StringCol(segchannel_itemsize),
        "lineage": tables.UInt16Col(),
        "progeny_id": tables.UInt64Col(),
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
    Re-label all non-zero pixels with the given label.
    """
    new_mask = (mask != 0) * label
    return new_mask


def relabel_mask_complete(mask, old_labels, new_labels):
    """
    For each label, replace it with a new label.
    Input a mask with labeled pixels.
    """
    # FIXME problem when different regions have the same label -- this can happen,
    # for example, when two cells are of a different lineage but the same
    # progeny id. Either disallow more than 1 lineage at a time,
    # or come up with a way to color all lineages without interference.
    # As is, each repeat of a progeny id will further increment previous iterations of that progeny id.
    new_mask = numpy.zeros(mask.shape, dtype=mask.dtype)
    for (o, n) in zip(old_labels, new_labels):
        new_mask += ((mask == o).astype(mask.dtype) * n).astype(mask.dtype)

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


def make_kymograph_from_centroids(table, is_top_row, trench_length):
    """
    Input:
        a. Dataframe slice corresponding to a single trench with a given
            - File
            - FOV
            - Trench number
            - Segmentation channel
        b. Whether the row of trenches is top or bottom
        c. The length of the trenches (to not count cell divisions for cells leaving tne trench)

    Output a "kymograph" stored as a list of lists:
        a. Each sub-list represents the cells at each subsequent time frame.
        b. Its contents are the cells, ordered from the feeding channel to the end of the trench:
           Two choices:
               - Always re-order (re-label) them by comparing centroids (safe)
               - or assume that cells were labeled in the correct order (risky)
            We choose to apply the sorting, using is_top_row as a guide.
        c. A "cell" is stored in said list as:
            - the connected component region's label (integer value).
            - the cell length
            - whether or not the cell overlaps the tench opening (the "edge")

    Also return the total number of cells.
    """
    kymograph = []
    total_cells = 0

    # Get unique frame numbers from the table, in sorted order
    # NOTE this code assumes no gaps (indexing by lists), and in sorted frame order
    # Gaps could arise, for example, for trenches at the edge of an FOV, which can fall out due to stage drift,
    # or if a frame went missing for any other reason.
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
        # NOTE: we use the "centroid_col" because the arrays are F-ordered.
        # A more typical C-ordered image would use the "centroid_row" to sort.
        cells_this_frame.sort_values(
            by=["centroid_col"], inplace=True, ascending=is_top_row
        )

        # Iterate in this new order and append the labels to a new list
        cell_order = []
        for c in cells_this_frame.itertuples():
            label = c.info_label
            length = c.axis_length_major

            # A bool to state whether a cell touches the trench opening
            # Useful for not erroneously calling cell divisions when a cell is exiting the trench
            if is_top_row:
                not_cell_edge = c.bounding_box_max_col != trench_length
            else:
                not_cell_edge = c.bounding_box_min_col != 0

            this_cell = {
                "label": label,
                "cell_length": length,
                "not_cell_edge": not_cell_edge,
            }

            cell_order.append(this_cell)
            total_cells += 1

        # And append this list to the kymograph
        kymograph.append(cell_order)

    return kymograph, total_cells


def run_kymograph_to_lineages(kymograph, total_cells, length_buffer):
    """
    Input
        - labeled kymograph for a given FOV / trench
        - total number of cells (used to initialize numpy arrays)
        - length buffer (in pixels) for calling cell divisions
          (allow small decreases in size to not be flagged as a division)

    Return a Pandas dataframe:
        - info_frame
        - info_label
        - lineage
        - progeny_id

    This dataframe can then be merged back with the pre-existing information (File, FOV, Z_level, etc.),
    thereby adding lineage & progeny_id columns (the lineage information).

    Each cell at frame 0 is said to be the progenitor of its own lineage.
    Every lineage gives rise to progeny.
    Progeny are assigned unique identifiers according to a binary tree system (n: 2n, 2n+1) (1: 2,3; 2: 4,5; 3: 6,7 etc.)

    This algorithm makes assumptions: see above.

    It works by greedily handling cells which are related to the deepest cell, given the above assumptions.
    Priority is handled by carefully managing the order in which cells are appended to & popped from the queue.
    We call each cell to be processed a "job," which conveniently is just a dict (hashmap), and the queue is a "stack" (a list).
    When a cell is popped from the stack, then its data are immediately added to result, which is a dict of numpy arrays.
    Then, its progeny are added to the queue.

    Each lineage is handled with its own queue.

    Finally, the results dict is converted to a Pandas dataframe.

    TODO allow the user to change other parameters?
    """
    # This is a dict of numpy arrays.
    # It will be converted to a dataframe & returned at the end.
    # NOTE if we are using groupby, do we have to specify all of these fields,
    # or are they automatically returned as the Index of the series?
    # We don't actually have access to them!
    result = {
        "info_frame": numpy.empty(total_cells, dtype=numpy.uint16),
        "info_label": numpy.empty(total_cells, dtype=numpy.uint16),
        "lineage": numpy.empty(total_cells, dtype=numpy.uint16),
        "progeny_id": numpy.empty(total_cells, dtype=numpy.uint64),
    }

    # In the zeroth time frame, define cells as lineage progenitors.
    # For each, run the lineage algorithm on the remaining kymograph time frame entries.
    lineage = 0
    current_cell_count = 0
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
            "not_cell_edge": cell["not_cell_edge"],
        }

        stack.append(job)

        # While it's not empty
        while stack:
            previous_cell = stack.pop()

            # Write the resulting cell to the dict, which will be
            # returned as a dataframe.
            result["info_frame"][current_cell_count] = previous_cell["frame"]
            result["info_label"][current_cell_count] = previous_cell["label"]
            result["lineage"][current_cell_count] = previous_cell["lineage"]
            result["progeny_id"][current_cell_count] = previous_cell["progeny_id"]
            current_cell_count += 1

            # Now look towards the next time frame
            this_progeny_id = previous_cell["progeny_id"]
            next_frame = previous_cell["frame"] + 1

            # Are there more frames to process?
            # Are there any cells in the next frame?
            if (next_frame < len(kymograph)) and (kymograph[next_frame]):
                # Load the next frame from the kymograph
                next_frame_group = kymograph[next_frame]

                # Get the top cell: contains a label and a major axis length
                top_cell = next_frame_group.pop(0)

                # Is the top cell shorter than the cell in the previous frame? -> There was a cell division
                # BUT if a cell borders the edge, then a change in size cannot be differentiated from
                # division or sliding out of the trench.
                # FIXME the check doesn't seem to work?
                if (
                    top_cell["cell_length"] + length_buffer
                    < previous_cell["cell_length"]
                ) and (job["not_cell_edge"]):

                    # If the cell is ~ 2x larger than that in the previous frame, then
                    # there might be a cell merging issue (segmentation problem).
                    # NOTE this didn't work. Clearly there are some mergers, but no warnings are printed.
                    if top_cell["cell_length"] > 2 * previous_cell["cell_length"]:
                        print("Warning: possible cell merger detected.")

                    # Is there also a next (sister) cell, from this division? Add it to the stack first.
                    # It's added first, because the stack is FILO.
                    if next_frame_group:
                        # Contains a label and a major axis length
                        next_cell = next_frame_group.pop(0)

                        job = {
                            "frame": next_frame,
                            "lineage": lineage,
                            "progeny_id": this_progeny_id * 2 + 1,
                            "label": next_cell["label"],
                            "cell_length": next_cell["cell_length"],
                            "not_cell_edge": cell["not_cell_edge"],
                        }

                        stack.append(job)

                    # Now add the top cell to the stack
                    job = {
                        "frame": next_frame,
                        "lineage": lineage,
                        "progeny_id": this_progeny_id * 2,
                        "label": top_cell["label"],
                        "cell_length": top_cell["cell_length"],
                        "not_cell_edge": cell["not_cell_edge"],
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
                        "not_cell_edge": cell["not_cell_edge"],
                    }

                    stack.append(job)

            # No cells: do nothing!

        # End while loop: stack is empty, can move on to the next lineage, with a new stack and any remaining cells left in the kymograph.
        lineage += 1

    # Finished, so convert the dict of arrays to a dataframe.
    result = pandas.DataFrame.from_dict(result)
    return result


# TODO it might be possible to walk through the dataframe using
# groupby and apply or map etc., something like:
def run_lineages(df, trench_length, length_buffer, row_number):
    # This is an entire column of row numbers. They are all the same.
    # Convert to a single number.
    # FIXME can I just use .loc  to get the zeroth element, without calling unique()?
    row_number = row_number.unique()[0]

    # Check the number of the row: even or odd, top or bottom
    # NOTE this assumes that there are only 2 rows!!
    # For other configurations, there would need to be a different way of specifying whether
    # the trench opening is at the bottom (as in a top row), or at the top (as in a bottom row).
    # Howvever, a boolean value could still be passed along.
    # The location of the trench opening is crucial for determinining the direction of cell division.
    is_top_row = (row_number % 2) == 0

    # Creates a "kymograph" (list of lists)
    (kymograph, total_cells) = make_kymograph_from_centroids(
        df, is_top_row, trench_length
    )

    # This returns a dataframe
    result = run_kymograph_to_lineages(kymograph, total_cells, length_buffer)
    return result


###
def main_lineages_function(infile, outfile, length_buffer, trench_length):
    """
    infile contains cell measurements table

    outfile contains a new table with the lineage information

    length_buffer for calling cell divisions, in pixels

    trench_length for flagging cells as touching the trench opening,
    in which case we disable calling cell divisions.
    """
    # TODO/FIXME infer trench length from the regions or somewhere else
    # trench_length = 152

    # Open the measurements table as a Pandas dataframe
    # NOTE this code assumes that the data were processed
    # after re-numbering the trenches to account for stage drift over time.
    # (see: renumber_trenches.py), producing a column called
    # corrected_trench_label.
    h5_in = tables.open_file(infile, "r")

    # NOTE This version of the code opts to load the entire
    # table into memory. Potential advatages include implicit
    # iteration with possible vectorization, and simpler
    # datatable merging & writing at the end.
    # Possible disadvantage is if the table doesn't fit into memory.
    # In that case, it would be necessary to iterate through chunks
    # using read_where() and to write dataframe chunks,
    # then to concatenate them all at the end.
    df = pandas.DataFrame(h5_in.get_node("/cell_measurements").read())
    h5_in.close()

    # Run lineage tracking & add the results as 2 new columns:
    # 1. lineage 2. progeny_id
    # Pass in the row number too so that we can check if it is top (even) or bottom (odd)
    # Return a new dataframe with a multi-index.
    grouping = [
        "info_file",
        "info_fov",
        "info_z_level",
        "info_row_number",
        "corrected_trench_label",
        "info_seg_channel",
    ]

    result = df.groupby(grouping).apply(
        lambda group: run_lineages(
            group, trench_length, length_buffer, group.info_row_number
        )
    )

    # Remove the multi-index.
    result = result.reset_index()
    # A new column called "level_6" is created from the within-group numbering.
    # We don't care about the numbering within groups, so drop it.
    result = result.drop(["level_6"], axis=1)

    # Merge into the original dataframe using file, fov, frame, z_level, row_number, corrected_trench_label, seg_channel, info_label
    on = [
        "info_file",
        "info_fov",
        "info_frame",
        "info_z_level",
        "info_row_number",
        "corrected_trench_label",
        "info_seg_channel",
        "info_label",
    ]

    df_merge = pandas.merge(df, result, on=on)
    # TODO is there a way for pandas to never convert bytes columns into string objects in the first place?
    df_merge["info_file"] = result["info_file"].astype(numpy.bytes_)
    df_merge["info_seg_channel"] = result["info_seg_channel"].astype(numpy.bytes_)
    write_dataframe_to_hdf5(df_merge, outfile, "cell_measurements")
