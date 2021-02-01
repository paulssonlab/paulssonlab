#!/usr/bin/env python
import os
import tables
import numpy
import napari
from dask import delayed
import dask.array
import pandas
from params import read_params_string, read_params_file
from napari_browse_hdf5 import (
    get_largest_extents_hdf5,
    metadata_attributes_equal,
    metadata_array_equal,
)
from lineages import relabel_mask_complete
import matplotlib.cm

"""
Browse kymographs from an HDF5 file using Napari.
"""


def dask_stack_loop(
    num_trenches,
    num_rows,
    file_nodes,
    fields_of_view,
    z_levels,
    tr_width,
    img_width,
    height,
    loading_function,
    dtype,
    load_func_args,
):
    """
    Create a "lazy" dask array to display information in Napari.
    This is a generic routine which should abstract across many different kinds of data.
    """
    lazy_func = delayed(loading_function)

    file_array = []
    for file in file_nodes:
        # FOV nodes
        fov_array = []
        for fov in fields_of_view:
            # Trench Rows
            row_array = []
            for row in range(num_rows):
                # Trenches
                tr_array = []
                for tr in range(num_trenches):
                    # z_levels
                    z_array = []
                    for z in z_levels:
                        arr = dask.array.from_delayed(
                            lazy_func(
                                tr,
                                row,
                                file,
                                fov,
                                z,
                                tr_width,
                                img_width,
                                height,
                                *load_func_args
                            ),
                            shape=(img_width, height),
                            dtype=dtype,
                        )
                        z_array.append(arr)

                    z_array_to_stack = dask.array.stack(z_array)
                    tr_array.append(z_array_to_stack)

                tr_array_to_stack = dask.array.stack(tr_array)
                row_array.append(tr_array_to_stack)

            row_array_to_stack = dask.array.stack(row_array)
            fov_array.append(row_array_to_stack)

        fov_array_to_stack = dask.array.stack(fov_array)
        file_array.append(fov_array_to_stack)

    megastack = dask.array.stack(file_array)

    return megastack


def query_regions_from_file(
    trench_number, row_number, file, fov, z_level, regions_file
):
    """
    Read all regions matching a given query
    Return a list of regions, and a list of trench numbers
    """
    h5file = tables.open_file(regions_file, "r")
    table = h5file.get_node("/regions_table")
    query = """(corrected_trench_label == trench_number) & \
               (info_row_number == row_number) & \
               (info_file == file) & \
               (info_fov == fov) & \
               (info_z_level == z_level)"""

    search = table.read_where(query)
    h5file.close()

    df = pandas.DataFrame(search)

    regions = {}
    for r in df.itertuples():
        # NOTE: have to add back the stage x, y positions,
        # since the original image is not "corrected" for stage drift.
        regions[r.info_frame] = [
            r.min_row + r.info_x,
            r.min_col - r.info_y,
            r.max_row + r.info_x,
            r.max_col - r.info_y,
        ]

    return regions


def kymo_masks(
    trench_number,
    row_number,
    file,
    fov,
    z_level,
    tr_width,
    img_width,
    height,
    sc,
    masks_file,
    frames,
    regions_file,
):
    """
    Returns a 2D numpy array with all of the kymograph timepoints written sequentially, from "left" to "right".
    """
    h5file = tables.open_file(regions_file, "r")
    table = h5file.get_node("/regions_table")
    query = """(info_file == file) & \
               (info_fov == fov) & \
               (info_z_level == z_level) & \
               (info_row_number == row_number) & \
               (corrected_trench_label == trench_number)"""

    search = table.read_where(query)
    df = pandas.DataFrame(search)

    kymograph_array = numpy.empty(shape=(img_width, height), dtype=numpy.uint16)
    h5file_masks = tables.open_file(masks_file, "r")

    # Trenches across frames (over time)
    for frame in frames:
        try:
            # FIXME remove the .iloc once I fix the table merging stuff?
            # Right now it returns an entry for each segmented cell.
            tr = df[(df["info_frame"] == frame)]["info_trench_number"].iloc[0]
            node_path = "/masks/{}/FOV_{}/Frame_{}/Z_{}/{}/{}/region_{}".format(
                file, fov, frame, z_level, sc, row_number, tr
            )
            node = h5file_masks.get_node(node_path)
            img_slice = node.read()
        except:
            img_slice = numpy.zeros(shape=(tr_width, height), dtype=numpy.uint16)

        kymograph_array[frame * tr_width : (frame + 1) * tr_width, 0:height] = img_slice

    h5file_masks.close()
    h5file.close()

    return kymograph_array.T


def kymo_mask_layer(
    num_trenches,
    num_rows,
    file_nodes,
    fields_of_view,
    frames,
    z_levels,
    tr_width,
    img_width,
    height,
    viewer,
    sc,
    masks_file,
    regions_file,
):
    """
    Generate kymographs of masks data.
    """
    for c in sc:
        mask_stack = dask_stack_loop(
            num_trenches,
            num_rows,
            file_nodes,
            fields_of_view,
            z_levels,
            tr_width,
            img_width,
            height,
            kymo_masks,
            numpy.uint16,
            [c, masks_file, frames, regions_file],
        )

        viewer.add_labels(mask_stack, name="Segmented Cells {}".format(sc))


def kymo_computed_image_layer(
    num_trenches,
    num_rows,
    file_nodes,
    fields_of_view,
    frames,
    z_levels,
    tr_width,
    img_width,
    height,
    viewer,
    masks_file,
    regions_file,
    data_table_file,
    viewer_params,
):
    """
    Generate computed images based on cell properties & masks.
    Input a column name from a data table, and use that to fill in the pixels of a mask,
    and present it as an intensity image.
    """
    for channel, categories in viewer_params.items():
        for category, items in categories.items():
            eval_string = items["eval_string"]
            params = items["viewer_params"]

            mask_stack = dask_stack_loop(
                num_trenches,
                num_rows,
                file_nodes,
                fields_of_view,
                z_levels,
                tr_width,
                img_width,
                height,
                kymo_computed_image,
                numpy.uint16,
                [
                    channel,
                    masks_file,
                    frames,
                    regions_file,
                    eval_string,
                    data_table_file,
                ],
            )

            viewer.add_image(mask_stack, **params)


def query_for_masks(file, fov, z, row_number, trench_number, sc, data_table_file):
    """
    Return a dataframe matching the query.
    Useful for getting matching cells in a cell measurements table.
    """
    h5file_data = tables.open_file(data_table_file, "r")
    data_table = h5file_data.get_node("/cell_measurements")
    # NOTE the speed of the queries might be sensitive to the order
    # in which the data are searched. If there are many trenches, but few files,
    # then it's much faster to search trenches first, then files.
    query = """(corrected_trench_label == trench_number) & \
               (info_row_number == row_number) & \
               (info_fov == fov) & \
               (info_file == file) & \
               (info_z_level == z) & \
               (info_seg_channel == sc)"""

    search = data_table.read_where(query)
    h5file_data.close()

    df = pandas.DataFrame(search)
    return df


def kymo_computed_image(
    trench_number,
    row_number,
    file,
    fov,
    z_level,
    tr_width,
    img_width,
    height,
    sc,
    masks_file,
    frames,
    regions_file,
    column_expression,
    data_table_file,
):
    """
    Returns a 2D numpy array with all of the kymograph timepoints written sequentially, from "left" to "right".
    This version is for filling in intensities by performing a query on the data table for specific quantities.
    """
    h5file = tables.open_file(regions_file, "r")
    table = h5file.get_node("/regions_table")
    query = """(corrected_trench_label == trench_number) & \
               (info_row_number == row_number) & \
               (info_file == file) & \
               (info_fov == fov) & \
               (info_z_level == z_level)"""

    search = table.read_where(query)
    df_regions = pandas.DataFrame(search)

    kymograph_array = numpy.empty(shape=(img_width, height), dtype=numpy.float32)
    h5file_masks = tables.open_file(masks_file, "r")

    df_this_row_trench = query_for_masks(
        file, fov, z_level, row_number, trench_number, sc, data_table_file
    )

    # Trenches across frames (over time)
    for frame in frames:
        try:
            # FIXME remove the .iloc once I fix the table merging stuff?
            # Right now it returns an entry for each segmented cell.
            tr = df_regions[(df_regions["info_frame"] == frame)][
                "info_trench_number"
            ].iloc[0]
            node_path = "/masks/{}/FOV_{}/Frame_{}/Z_{}/{}/{}/region_{}".format(
                file, fov, frame, z_level, sc, row_number, tr
            )
            node = h5file_masks.get_node(node_path)
            mask = node.read()

            df_this_frame = df_this_row_trench[
                df_this_row_trench["info_frame"] == frame
            ]
            computed_image = numpy.zeros(mask.shape, dtype=numpy.float32)
            # Speed up the compuation by only masking the bbox area
            for cell in df_this_frame.itertuples():
                mask_area = mask[
                    cell.bounding_box_min_row : cell.bounding_box_max_row,
                    cell.bounding_box_min_col : cell.bounding_box_max_col,
                ]

                computed_image[
                    cell.bounding_box_min_row : cell.bounding_box_max_row,
                    cell.bounding_box_min_col : cell.bounding_box_max_col,
                ] = (mask_area == cell.info_label) * eval(column_expression)

        except:
            computed_image = numpy.zeros(shape=(tr_width, height), dtype=numpy.float32)

        kymograph_array[
            frame * tr_width : (frame + 1) * tr_width, 0:height
        ] = computed_image

    h5file_masks.close()
    h5file.close()

    return kymograph_array.T


def make_tab20_mapping(num_colors=65535):
    """
    Use the "tab20" colormap to color pairs of cells,
    and add "black" (transparent) for 0, and white for the first cell.

    tab20 is good because the colors come in pairs with similar hues,
    whereas the default mapping in Napari is more "randomized," which
    makes it hard to visually match pairs of cells.

    NOTE/FIXME: how high do we want the labels to be? Use 65535 for now,
    but since the labels increase as powers of 2, this might not be enough
    for very long timelapses.
    Could search the lineages file & find the maximum possible label?
    """
    tab20 = matplotlib.cm.get_cmap("tab20")

    mapping = {}
    # Transparent for background
    mapping[0] = (0.0, 0.0, 0.0)
    # White for the first cell
    mapping[1] = (1.0, 1.0, 1.0)
    for i in range(num_colors):
        # The tab20 map has 20 colors, so cycle them every 20 using modulo
        mapping[i + 2] = tab20((i + 2) % 20)

    return mapping


def kymo_lineages_layer(
    num_trenches,
    num_rows,
    file_nodes,
    fields_of_view,
    frames,
    z_levels,
    tr_width,
    img_width,
    height,
    viewer,
    sc,
    masks_file,
    lineages_file,
):
    """
    Generate kymographs of masks data and re-label the cells based on a computed lineage.
    """
    for c in sc:
        c = c.decode("utf-8")
        mask_stack = dask_stack_loop(
            num_trenches,
            num_rows,
            file_nodes,
            fields_of_view,
            z_levels,
            tr_width,
            img_width,
            height,
            kymo_lineage_masks,
            numpy.uint16,
            [c, frames, masks_file, lineages_file],
        )

        mapping = make_tab20_mapping()
        viewer.add_labels(mask_stack, name="Lineages ({})".format(c), color=mapping)


def kymo_lineage_masks(
    trench_number,
    row_number,
    file,
    fov,
    z_level,
    tr_width,
    img_width,
    height,
    sc,
    frames,
    masks_file,
    lineages_file,
):
    """
    Returns a 2D numpy array with all of the kymograph timepoints written sequentially, from "left" to "right".
    In addition, re-label the cells based on information from a computed lineage.
    """
    h5file_lineages = tables.open_file(lineages_file, "r")
    table_lineages = h5file_lineages.get_node("/cell_measurements")

    query = """(corrected_trench_label == trench_number) & \
               (info_row_number == row_number) & \
               (info_file == file) & \
               (info_fov == fov) & \
               (info_z_level == z_level) & \
               (info_seg_channel == sc)"""

    cells_of_interest = table_lineages.read_where(query)
    df_cells = pandas.DataFrame(cells_of_interest)

    kymograph_array = numpy.empty(shape=(img_width, height), dtype=numpy.uint16)
    h5file_masks = tables.open_file(masks_file, "r")

    # Trenches across frames (over time)
    for frame in frames:
        new_mask = numpy.zeros(shape=(tr_width, height), dtype=numpy.uint16)
        try:
            # FIXME remove the .iloc once I fix the table merging stuff?
            # Right now it returns an entry for each segmented cell.
            this_frame = df_cells[df_cells["info_frame"] == frame]
            tr = this_frame["info_trench_number"].iloc[0]
            node_path = "/masks/{}/FOV_{}/Frame_{}/Z_{}/{}/{}/region_{}".format(
                file, fov, frame, z_level, sc, row_number, tr
            )
            node = h5file_masks.get_node(node_path)
            mask = node.read()

            # Re-label the labeled mask slice
            all_lineages = this_frame["lineage"].unique()
            for lineage in all_lineages:
                # Re-label this lineage according to the progeny id
                this_lineage = this_frame[this_frame["lineage"] == lineage]
                old_labels = this_lineage["info_label"].values

                # Run a series of logical ORs across all labels to get all matching pixels
                # This "mask" only contains pixels from this lineage.
                this_lineage_mask = numpy.zeros(new_mask.shape, dtype=numpy.bool)
                for o in old_labels:
                    this_lineage_mask |= mask == o

                # Convert back to the original dtype so that we can do sums & multiplications
                this_lineage_mask = this_lineage_mask.astype(new_mask.dtype)

                # Apply the lineage mask to the original labels mask
                labeled_pixels_this_lineage = mask * this_lineage_mask

                # Now get the progeny_id as a mask too:
                # Replace the labels in the lineage pixels with their progeny_id
                new_labels = this_lineage["progeny_id"].values
                labeled_pixels_progeny = relabel_mask_complete(
                    labeled_pixels_this_lineage, old_labels, new_labels
                )
                new_mask += labeled_pixels_progeny

        except Exception as e:
            print(str(e))

        kymograph_array[frame * tr_width : (frame + 1) * tr_width, 0:height] = new_mask

    h5file_masks.close()
    h5file_lineages.close()

    return kymograph_array.T


def kymo_img(
    trench_number,
    row_number,
    file,
    fov,
    z_level,
    tr_width,
    img_width,
    height,
    channel,
    regions_file,
    images_file,
    frames,
):
    """
    Returns a 2D numpy array with all of the kymograph timepoints written sequentially, from "left" to "right".
    """
    matching_frames = query_regions_from_file(
        trench_number, row_number, file, fov, z_level, regions_file
    )
    kymograph_array = numpy.empty(shape=(img_width, height), dtype=numpy.float32)
    h5file_images = tables.open_file(images_file, "r")

    # Trenches across frames (over time)
    for fr in frames:
        # NOTE: Can this ever fail?
        try:
            tr = matching_frames[fr]
            node_path = "/Images/{}/FOV_{}/Frame_{}/Z_{}/{}".format(
                file, fov, fr, z_level, channel
            )
            node = h5file_images.get_node(node_path)
            img_slice = node[tr[0] : tr[2], tr[1] : tr[3]].astype(numpy.float32)
            # TODO: requires also taking a slice out of the corrections array
        except:
            img_slice = numpy.zeros(shape=(tr_width, height), dtype=numpy.float32)

        kymograph_array[fr * tr_width : (fr + 1) * tr_width, 0:height] = img_slice

    h5file_images.close()

    return kymograph_array.T


def kymo_img_layers(
    num_trenches,
    num_rows,
    file_nodes,
    fields_of_view,
    frames,
    z_levels,
    viewer,
    tr_width,
    img_width,
    height,
    channels,
    layer_params,
    corrections_file,
    regions_file,
    images_file,
):
    """
    Generate kymographs of image data. TODO: Optionally apply flatfield & other corrections.
    """
    for c in channels:
        image_stack = dask_stack_loop(
            num_trenches,
            num_rows,
            file_nodes,
            fields_of_view,
            z_levels,
            tr_width,
            img_width,
            height,
            kymo_img,
            numpy.float32,
            [c, regions_file, images_file, frames],
        )  # TODO: corrections, corrections_file, c])

        viewer.add_image(image_stack, **layer_params[c])


def get_reg_dims(regions_file):
    """
    Use the regions h5file table to infer trench dimensions,
    max number of rows of trenches, and max number of trenches per FOV.
    """
    h5file_reg = tables.open_file(regions_file, "r")
    table = h5file_reg.get_node("/regions_table")
    df_reg = pandas.DataFrame(table.read())
    h5file_reg.close()

    # Assume all trenches have the same dimensions,
    # and use the zeroth trench to determine what they are.
    first_element = df_reg.iloc[0]
    height = first_element.max_col - first_element.min_col
    tr_width = first_element.max_row - first_element.min_row

    num_rows = df_reg["info_row_number"].max() + 1

    # This requires running trench drift correction
    num_trenches = df_reg["corrected_trench_label"].max()

    return (tr_width, height, num_rows, num_trenches)


def main_kymograph_browser_function(
    images_file,
    masks_file,
    regions_file,
    napari_settings_file,
    corrections_file,
    lineages_file,
    viewer_params_file,
    data_table_file,
):
    """
    Main function, invoked from the command-line.
    """
    with napari.gui_qt():
        viewer = napari.Viewer()

        # Load the image layer params for napari
        # TODO params for masks layers? regions layers?
        layer_params = read_params_file(napari_settings_file)

        # File with images & metadata
        h5file = tables.open_file(images_file, "r")

        # Get the largest extents across all the nd2 files within the h5 file,
        # so that the dask array is made with the proper min, max dimensions for each dimension.
        attributes = ["fields_of_view", "frames", "z_levels"]
        extents = {}
        for a in attributes:
            (smallest, largest) = get_largest_extents_hdf5(h5file, a)
            if smallest == largest:
                extents[a] = [smallest]
            else:
                extents[a] = [i for i in range(smallest, largest + 1)]

        # Define the channels
        channels = metadata_array_equal(h5file, "channels")
        channels = [c.decode("utf-8") for c in channels]

        # Iterate the H5 file images & lazily load them into a dask array
        file_nodes = [x._v_name for x in h5file.list_nodes("/Images")]

        h5file.close()

        # Open the regions table, get the first region, and
        # recalculate based on min_row, max_row; min_col, max_col.
        (tr_width, height, num_rows, num_trenches) = get_reg_dims(regions_file)

        # This is the max value of info_frame. Or unique.
        num_frames = len(extents["frames"])

        # Use the above to calculate this. Used to create empty canvases.
        img_width = num_frames * tr_width

        # Image layers
        kymo_img_layers(
            num_trenches=num_trenches,
            num_rows=num_rows,
            file_nodes=file_nodes,
            fields_of_view=extents["fields_of_view"],
            frames=extents["frames"],
            z_levels=extents["z_levels"],
            viewer=viewer,
            tr_width=tr_width,
            img_width=img_width,
            height=height,
            channels=channels,
            layer_params=layer_params,
            corrections_file=corrections_file,
            regions_file=regions_file,
            images_file=images_file,
        )

        # Masks layers
        if masks_file:
            masks_file_path = os.path.join(masks_file, "MASKS/masks.h5")

            # Load the segmentation channels from the masks h5file
            h5file_masks = tables.open_file(masks_file_path, "r")

            seg_params_node = h5file_masks.get_node("/Parameters", "seg_params.yaml")
            seg_params_str = seg_params_node.read().tobytes().decode("utf-8")
            seg_params_dict = read_params_string(seg_params_str)
            seg_channels = [k for k in seg_params_dict.keys()]

            h5file_masks.close()

            kymo_mask_layer(
                num_trenches=num_trenches,
                num_rows=num_rows,
                file_nodes=file_nodes,
                fields_of_view=extents["fields_of_view"],
                frames=extents["frames"],
                z_levels=extents["z_levels"],
                tr_width=tr_width,
                img_width=img_width,
                height=height,
                viewer=viewer,
                sc=seg_channels,
                masks_file=masks_file_path,
                regions_file=regions_file,
            )

        # Lineages layers
        if lineages_file:
            # Load the segmentation channels from the measurements & lineages table
            h5file_lineages = tables.open_file(lineages_file, "r")
            table = h5file_lineages.get_node("/cell_measurements")
            df = pandas.DataFrame(table.read())
            seg_channels = df["info_seg_channel"].unique()
            # TODO what about converting to / from bytes & utf-8?
            h5file_lineages.close()
            masks_file_path = os.path.join(masks_file, "MASKS/masks.h5")

            kymo_lineages_layer(
                num_trenches=num_trenches,
                num_rows=num_rows,
                file_nodes=file_nodes,
                fields_of_view=extents["fields_of_view"],
                frames=extents["frames"],
                z_levels=extents["z_levels"],
                tr_width=tr_width,
                img_width=img_width,
                height=height,
                viewer=viewer,
                sc=seg_channels,
                masks_file=masks_file_path,
                lineages_file=lineages_file,
            )

        # Input a set of images, but also run a calculation on them before displaying.
        # Store the result of the calculation within the regions defined by cell masks.
        if viewer_params_file and data_table_file and masks_file:
            viewer_params = read_params_file(viewer_params_file)
            masks_file_path = os.path.join(masks_file, "MASKS/masks.h5")

            kymo_computed_image_layer(
                num_trenches=num_trenches,
                num_rows=num_rows,
                file_nodes=file_nodes,
                fields_of_view=extents["fields_of_view"],
                frames=extents["frames"],
                z_levels=extents["z_levels"],
                tr_width=tr_width,
                img_width=img_width,
                height=height,
                viewer=viewer,
                masks_file=masks_file_path,
                regions_file=regions_file,
                data_table_file=data_table_file,
                viewer_params=viewer_params,
            )
