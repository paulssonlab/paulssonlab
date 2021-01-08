#!/usr/bin/env python

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
from lineages import relabel_mask_complete, relabel_mask
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
    # df = pandas.read_hdf(regions_file, "/regions_table")
    h5file = tables.open_file(regions_file, "r")
    table = h5file.get_node("/regions_table")
    query = """(info_file == file) & \
               (info_fov == fov) & \
               (info_z_level == z_level) & \
               (info_row_number == row_number) & \
               (corrected_trench_label == trench_number)"""

    search = table.read_where(query)
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

    h5file.close()

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
    sc,
    masks_file,
    regions_file,
    computed_image_function,
    data_table_file,
    data_table_name,
):
    """
    Generate computed images based on cell properties & masks.
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
            kymo_computed_image,
            numpy.uint16,
            [
                c,
                masks_file,
                frames,
                regions_file,
                computed_image_function,
                data_table_file,
                data_table_name,
            ],
        )

        viewer.add_image(mask_stack, name="Computed Image {}".format(c))


# TODO: select a function from a list of available functions,
# specified on the command line.
def run_computed_image_function(mask, df):
    # Example function: for every unique cell label, replace its
    # contents with the mean of its mVenus fluorescence intensities.
    # Calculate mVenus / pixel
    # Make a copy if we want to add a new column.
    df = df.copy()
    df["mean_mVenus"] = df["total_intensity_YFP"] / df["geometry_Area"]
    key = "mean_mVenus"

    # Alternatively: just display the Area without calculating anything
    # key = "geometry_Area"

    new_img = numpy.zeros(mask.shape, dtype=numpy.float32)
    for cell in df.itertuples():
        # TODO Can this be sped up by doing pixel_value -> f(pixel_value),
        # where f() is input cell label, output a const number?
        # Avoids looping. But it's actually pretty fast now.
        new_img += (mask == cell.info_label) * getattr(cell, key)
        # DEBUG
        # new_img = mask.astype(numpy.float32)

    return new_img


def query_for_masks(
    file, fov, z, row_number, trench_number, sc, data_table_file, data_table_name
):
    """
    Return a dataframe matching the query.
    Useful for getting matching cells in a cell measurements table.
    """
    h5file_data = tables.open_file(data_table_file, "r")
    data_table = h5file_data.get_node("/", data_table_name)
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
    computed_image_function,
    data_table_file,
    data_table_name,
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
    df_regions = pandas.DataFrame(search)

    kymograph_array = numpy.empty(shape=(img_width, height), dtype=numpy.float32)
    h5file_masks = tables.open_file(masks_file, "r")

    df_this_row_trench = query_for_masks(
        file,
        fov,
        z_level,
        row_number,
        trench_number,
        sc,
        data_table_file,
        data_table_name,
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
            computed_image = computed_image_function(mask, df_this_row_trench)

        except:
            computed_image = numpy.zeros(shape=(tr_width, height), dtype=numpy.float32)

        kymograph_array[
            frame * tr_width : (frame + 1) * tr_width, 0:height
        ] = computed_image

    h5file_masks.close()
    h5file.close()

    return kymograph_array.T


def make_tab20_mapping(num_colors=4096):
    """
    Use the "tab20" colormap to color pairs of cells,
    and add "black" (transparent) for 0, and white for the first cell.
    NOTE/FIXME: how high do we want the labels to be? Use 4096 for now,
    but since the labels increase as powers of 2, this might not be enough
    for very long timelapses.
    tab20 is good because the colors come in pairs with similar hues,
    whereas the default mapping in Napari is more "randomized," which
    makes it hard to visually match pairs of cells.
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
    # NOTE thus far lineages weren't computed using z-levels.
    """
    h5file_lineages = tables.open_file(lineages_file, "r")
    table_lineages = h5file_lineages.get_node("/cells_and_lineages")

    query = """(info_file == file) & \
               (info_fov == fov) & \
               (info_z_level == z_level) & \
               (info_seg_channel == sc) & \
               (info_row_number == row_number) & \
               (corrected_trench_label == trench_number)"""

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
            tr = this_frame["info_trench_number"].iloc[
                0
            ]  # FIXME sometimes getting out-of-bounds??
            node_path = "/masks/{}/FOV_{}/Frame_{}/Z_{}/{}/{}/region_{}".format(
                file, fov, frame, z_level, sc, row_number, tr
            )
            node = h5file_masks.get_node(node_path)
            mask = node.read()

            # Re-label the labeled mask slice
            all_lineages = this_frame["lineage"].unique()
            all_lineages = [0]
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

                # TODO purge OLD CODE
                ## Filter out all other lineages from the mask
                ## If it matches this lineage, replace the lineage with a 1
                ## Otherwise, replace it with a 0.
                # mask_this_lineage = relabel_mask(mask, lineage, 1, 0).astype(numpy.bool)

                ## Apply this filter to the original mask -> only keep pixels from this lineage,
                ## but preserve their original labels.
                # mask_this_lineage = mask * mask_this_lineage

                ## Re-label this lineage according to the progeny id
                # this_lineage = this_frame[this_frame["lineage"] == lineage]
                # old_labels = this_lineage["info_label"].values
                # new_labels = this_lineage["progeny_id"].values
                # new_mask += relabel_mask_complete(mask_this_lineage, old_labels, new_labels)
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


def main_kymograph_browser_function(
    images_file,
    masks_file,
    regions_file,
    napari_settings_file,
    corrections_file,
    lineages_file,
    computed_image_function,
    data_table_file,
    data_table_name,
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

        # Get the height & width, while checking that they are identical across all nd2 files
        height = metadata_attributes_equal(h5file, "height")
        width = metadata_attributes_equal(h5file, "width")

        # Define the channels
        channels = metadata_array_equal(h5file, "channels")
        channels = [c.decode("utf-8") for c in channels]

        # Iterate the H5 file images & lazily load them into a dask array
        file_nodes = [x._v_name for x in h5file.list_nodes("/Images")]

        h5file.close()

        # FIXME Need to open the merged cell & region info table, with
        # corrected trench bboxes.

        # Or, could open the regions table, get the first region, and
        # recalculate based on min_row, max_row; min_col, max_col.
        tr_width = 15
        height = 153

        # What to do if this value changes from one FOV to another?
        # Need to set this to the max of all possible (or unique?), and then
        # fill in empty as necessary. Can extract this
        num_rows = 2

        # This needs to be the max value of the drift-corrected,
        # re-labeled trenches.
        num_trenches = 96

        # This is the max value of info_frame. Or unique.
        num_frames = len(extents["frames"])

        # Use the above to calculate this. Used to create empty canvases.
        img_width = num_frames * tr_width

        # Image layers
        kymo_img_layers(
            num_trenches,
            num_rows,
            file_nodes,
            extents["fields_of_view"],
            extents["frames"],
            extents["z_levels"],
            viewer,
            tr_width,
            img_width,
            height,
            channels,
            layer_params,
            corrections_file,
            regions_file,
            images_file,
        )

        # Masks layers
        if masks_file:
            h5file_masks = tables.open_file(masks_file, "r")

            # Load the segmentation channels from the masks h5file
            seg_params_node = h5file_masks.get_node("/Parameters", "seg_params.yaml")
            seg_params_str = seg_params_node.read().tobytes().decode("utf-8")
            seg_params_dict = read_params_string(seg_params_str)
            seg_channels = [k for k in seg_params_dict.keys()]

            h5file_masks.close()

            kymo_mask_layer(
                num_trenches,
                num_rows,
                file_nodes,
                extents["fields_of_view"],
                extents["frames"],
                extents["z_levels"],
                tr_width,
                img_width,
                height,
                viewer,
                seg_channels,
                masks_file,
                regions_file,
            )

        # Lineages layers
        if lineages_file:
            # Load the segmentation channels from the measurements & lineages table
            h5file_lineages = tables.open_file(lineages_file, "r")
            table = h5file_lineages.get_node("/cells_and_lineages")
            df = pandas.DataFrame(table.read())
            seg_channels = df["info_seg_channel"].unique()
            # TODO what about converting to / from bytes & utf-8?
            h5file_lineages.close()

            kymo_lineages_layer(
                num_trenches,
                num_rows,
                file_nodes,
                extents["fields_of_view"],
                extents["frames"],
                extents["z_levels"],
                tr_width,
                img_width,
                height,
                viewer,
                seg_channels,
                masks_file,
                lineages_file,
            )

        # FIXME TEMP
        computed_image_function = run_computed_image_function
        # Would it make more sense to pre-compute the columns? Basically
        # make it a column display, and just pass in the name of the column.
        # Input a set of images, but also run a calculation on them before displaying.
        if computed_image_function and data_table_file and data_table_name:
            # TODO: specify which channels to perform the calculation
            # Pass in a list... best way? comma-separated? YAML?
            # These must match segmentation names, NOT microscope channel names.
            seg_channels = ["mKate2"]

            kymo_computed_image_layer(
                num_trenches,
                num_rows,
                file_nodes,
                extents["fields_of_view"],
                extents["frames"],
                extents["z_levels"],
                tr_width,
                img_width,
                height,
                viewer,
                seg_channels,
                masks_file,
                regions_file,
                computed_image_function,
                data_table_file,
                data_table_name,
            )
