import tables
import numpy
import napari
from dask import delayed
import dask.array
import pandas
from params import read_params_string, read_params_file
import os

# For registration
from imwarp import imwarp
from imref import imref2d

"""
NOTE for now, assumes image dtype is float32

TODO clean up the metadata functions, and their names, then stick them in the separate metadata Python file
Might also need to rework the similar functions used for nd2 browsing?

TODO 3D (z-stacked) layers need to add a scaling, as the Z dimension mightn't be the same as the X/Y.
X/Y pixel_microns is easy to grab from the metadata.
Z can be found in "/Metadata/File_beads_010/raw_metadata/image_metadata/SLxExperiment/uLoopPars/dZStep",
although I don't know how reliable this is.
The safest might be to always set the scale for X & Y using pixel microns,
and then try to load that value & use it to set the scale for Z.
"""


def dask_stack_loop(
    file_nodes,
    fields_of_view,
    frames,
    z_levels,
    width,
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
            # Frame nodes
            frame_array = []
            for frame in frames:
                z_array = []
                # z_levels
                for z in z_levels:
                    arr = dask.array.from_delayed(
                        lazy_func(file, fov, frame, z, width, height, *load_func_args),
                        shape=(width, height),
                        dtype=dtype,
                    )
                    z_array.append(arr)

                z_array_to_stack = dask.array.stack(z_array)
                frame_array.append(z_array_to_stack)

            frame_array_to_stack = dask.array.stack(frame_array)
            fov_array.append(frame_array_to_stack)

        fov_array_to_stack = dask.array.stack(fov_array)
        file_array.append(fov_array_to_stack)

    megastack = dask.array.stack(file_array)

    return megastack


### Regions (e.g. Trench Rows, Trenches) Layer


def add_regions_layer(
    regions_file, file_nodes, fields_of_view, frames, z_levels, viewer, width, height
):
    """
    There are 2 kinds of regions:
    1. Trench rows
    2. Trenches

    The individual elements (rows, trenches) must be pasted into a canvas which spans the original image.

    Create one label layer for trench rows, and another label layer for trenches.
    """
    # 1. Trench rows
    trench_rows_stack = dask_stack_loop(
        file_nodes,
        fields_of_view,
        frames,
        z_levels,
        width,
        height,
        regions_arr_to_canvas,
        numpy.uint16,
        [regions_file],
    )

    viewer.add_labels(trench_rows_stack, name="Trench Rows")

    # 2. Trenches within each row
    trenches_stack = dask_stack_loop(
        file_nodes,
        fields_of_view,
        frames,
        z_levels,
        width,
        height,
        trenches_arr_to_canvas,
        numpy.uint16,
        [regions_file],
    )

    viewer.add_labels(trenches_stack, name="Trenches")


def query_regions_from_file(file, fov, frame, z_level, regions_file):
    """
    Read all regions matching a given query
    Return a list of regions, and a list of trench numbers
    """
    h5file_regions = tables.open_file(regions_file, "r")
    regions_table = h5file_regions.get_node("/", "trench_coords")

    query = "(info_file == file) & \
              (info_fov == fov) & \
              (info_frame == frame) & \
              (info_z_level == z_level)"

    search = regions_table.read_where(query)
    df = pandas.DataFrame(search)

    h5file_regions.close()

    regions = []
    trench_numbers = []
    for r in df.itertuples():
        regions.append([r.min_row, r.min_col, r.max_row, r.max_col])
        trench_numbers.append(r.info_trench_number)

    return (regions, trench_numbers)


def trenches_arr_to_canvas(file, fov, frame, z_level, width, height, regions_file):
    """
    Read in bounding boxes for all the trenches, and write them to a canvas.
    """
    (regions, trench_numbers) = query_regions_from_file(
        file, fov, frame, z_level, regions_file
    )

    canvas = numpy.zeros(shape=(width, height), dtype=numpy.uint16)

    for r, tr in zip(regions, trench_numbers):
        # +1 because zero is background
        canvas[r[0] : r[2], r[1] : r[3]] = tr + 1

    return canvas.T


def regions_to_list(file, fov, frame, z_level, regions_file):
    """
    Return a list of region min_row, min_col, max_row, max_col matching the given condition.
    """
    # TODO: trench row number?
    h5file_regions = tables.open_file(regions_file, "r")
    regions_table = h5file_regions.get_node("/", "row_coords")

    condition = "(info_file == file) & \
                 (info_fov == fov) & \
                 (info_frame == frame) & \
                 (info_z_level == z_level)"

    search = regions_table.read_where(condition)
    df = pandas.DataFrame(search)
    h5file_regions.close()

    regions = []
    for r in df.itertuples():
        regions.append([r.min_row, r.min_col, r.max_row, r.max_col])

    return regions


def regions_arr_to_canvas(file, fov, frame, z_level, width, height, regions_file):
    """
    Input a numpy array with many rows, each containing the rectangular coordinates of regions.

    NOTE now, the input is a table containing relevant data,
    and then the table is queried to create the regions arr.

    Pass in info from which to read corresponding trench coords from the global trench coords table.

    Return a numpy array with numbered pixel positions for each region, and zeros for background.
    NOTE: what happens if the regions overlap?
    """
    regions = regions_to_list(file, fov, frame, z_level, regions_file)

    canvas = numpy.zeros(shape=(width, height), dtype=numpy.uint16)
    for i, r in enumerate(regions):
        canvas[r[0] : r[2], r[1] : r[3]] = i + 1

    return canvas.T


### Masks (e.g. cells) Layer


def add_masks_layer(
    masks_file,
    regions_file,
    file_nodes,
    fields_of_view,
    frames,
    z_levels,
    seg_channels,
    width,
    height,
    viewer,
):
    """
    Because each region gets its own masks, they need to be merged!
    This means that it is necessary to specify which regions were used to do so during the segmentation step.

    TODO figure out how to require regions? (or, the same region will just be over-written, and the user will figure it out?)

    As with the trenches, these regions the masks within these regions need to be combined into a single, large image.

    TODO pass in information for the labels layer params? name etc.
    """
    # If there were regions, then the masks need to be agglomerated into a larger image
    # Compile a canvas with all the masks from all the segmented areas
    # FIXME Would it make sense to use regions_arr_to_canvas() instead??
    # I think what this code is doing is recursing through all of the arrays on disk.
    # However, the point of the table was to replace all that with queries instead.
    # We can just iterate through the unique keys in the table,
    # and query them for the coordinates.
    lazy_masks_to_canvas = delayed(masks_to_canvas)

    # FIXME how to determine the number of rows? should each row set get its own labels layer?
    # that would mean sc x rows total number of layers.
    #
    # NOTE that if the regions within a set overlap, then they will overwrite each other
    # when it comes time to merging them together in the same canvas.
    if regions_file:
        for sc in seg_channels:
            masks_stack = dask_stack_loop(
                file_nodes,
                fields_of_view,
                frames,
                z_levels,
                width,
                height,
                masks_to_canvas,
                numpy.uint16,
                [sc, regions_file, masks_file],
            )

            viewer.add_labels(masks_stack, name="Segmented Cells {}".format(sc))

    # No regions file
    else:
        for sc in seg_channels:
            masks_stack = dask_stack_loop(
                file_nodes,
                fields_of_view,
                frames,
                z_levels,
                width,
                height,
                whole_labeled_canvas,
                numpy.uint16,
                [sc, masks_file],
            )

            viewer.add_labels(masks_stack, name="Segmented Cells {}".format(sc))

    # Done!


### Misc functions

# def load_correction_values(in_file):
# """
# Load values from an HDF5 file into a dict. Useful for flatfield corrections & camera biases.
# """
# if in_file:
# dict = {}

# h5 = tables.open_file(in_file, "r")
# for node in h5.iter_nodes("/"):
# key = node._v_name
# dict[key] = node.read()
# h5.close()

# else:
# dict = None

# return dict


def whole_labeled_canvas(file, fov, frame, z, width, height, sc, masks_file):
    """
    Return a labeled canvas, with no regions. This is done by just loading region_0.

    FIXME Passing in width & height, even though they aren't needed, because it conforms
    to the generic shape shared by all the other functions. There must be a better way to do this!
    Maybe by using named arguments, or a dict of arguments?
    """
    h5file_masks = tables.open_file(masks_file, "r")
    # If no regions, then default to the zeroth region set number and region
    path = "/masks/{}/FOV_{}/Frame_{}/Z_{}/{}/0/region_0".format(
        file, fov, frame, z, sc
    )
    mask = h5file_masks.get_node(path).read()
    h5file_masks.close()

    return mask.T


def regions_to_dict(regions_file, file, fov, frame, z_level):
    """
    Returns a dictionary mapping a region number to the a tuple of min_row, min_col, max_row, max_col.
    """
    h5file_regions = tables.open_file(regions_file, "r")
    table_regions = h5file_regions.get_node("/trench_coords")
    condition = "(info_file == file) & (info_fov == fov) & (info_frame == frame) & (info_z_level == z_level)"
    search_regions = table_regions.read_where(condition)
    h5file_regions.close()
    df_regions = pandas.DataFrame(search_regions)

    # Create a dict: input row number & trench number,
    # output bbox coords of the corresponding trench.
    regions = {}
    for r in df_regions.itertuples():
        key = "{}/region_{}".format(r.info_row_number, r.info_trench_number)
        regions[key] = (r.min_row, r.min_col, r.max_row, r.max_col)

    return regions


def masks_to_canvas(file, fov, frame, z, width, height, sc, regions_file, masks_file):
    """
    TODO update documentation how this function works
    Input a numpy array with many rows, each containing the rectangular coordinates of regions.

    Input a dict with keys structured as "rownum/region_regionnum", each containing the rectangular coordinates of regions.

    Return a numpy array with labeled pixel positions across all the masks, and zeros for background.
    NOTE: what happens if the regions overlap?
    """
    regions = regions_to_dict(regions_file, file, fov, frame, z)

    h5file_masks = tables.open_file(masks_file, "r")

    # Go through series of de-references b/c pytables doesn't allow iterating
    # over nodes if some of the nodes are symlinks.
    path = "/masks/{}".format(file)
    node_1 = h5file_masks.get_node(path)()
    node_2 = h5file_masks.get_node(node_1, "/FOV_{}".format(fov))()
    node_3 = h5file_masks.get_node(node_2, "/Frame_{}".format(frame))()
    node_4 = h5file_masks.get_node(node_3, "/Z_{}/{}".format(z, sc))

    canvas = numpy.zeros(shape=(width, height), dtype=numpy.uint16)
    for n in h5file_masks.walk_nodes(node_4, classname="CArray"):
        # Parse n's parent & its name to get keys
        row = n._v_parent._v_name
        trench = n.name.split("_")[1]
        key = "{}/region_{}".format(row, trench)
        r = regions[key]
        canvas[r[0] : r[2], r[1] : r[3]] = n.read()

    h5file_masks.close()

    return canvas.T


def query_for_masks(file, fov, frame, z, sc, data_table_file):
    """
    Return a dataframe matching the query.
    Useful for getting matching cells in a cell measurements table.
    """
    h5file_data = tables.open_file(data_table_file, "r")
    data_table = h5file_data.get_node("/cell_measurements")
    # NOTE the speed of the queries might be sensitive to the order
    # in which the data are searched. If there are many trenches, but few files,
    # then it's much faster to search trenches first, then files.
    query = """(info_frame == frame) & \
               (info_fov == fov) & \
               (info_file == file) & \
               (info_z_level == z) & \
               (info_seg_channel == sc)"""

    search = data_table.read_where(query)
    h5file_data.close()

    df = pandas.DataFrame(search)
    return df


def computed_masks_to_canvas(
    file,
    fov,
    frame,
    z,
    width,
    height,
    sc,
    regions_file,
    masks_file,
    column_expression,
    data_table_file,
):
    """
    TODO update documentation how this function works

    TODO can we speed up the compuation by only masking the bbox area?
    """
    regions = regions_to_dict(regions_file, file, fov, frame, z)

    df = query_for_masks(file, fov, frame, z, sc, data_table_file)

    # Go through series of de-references b/c pytables doesn't allow iterating
    # over nodes if some of the nodes are symlinks.
    h5file_masks = tables.open_file(masks_file, "r")
    path = "/masks/{}".format(file)
    node_1 = h5file_masks.get_node(path)()
    node_2 = h5file_masks.get_node(node_1, "/FOV_{}".format(fov))()
    node_3 = h5file_masks.get_node(node_2, "/Frame_{}".format(frame))()
    node_4 = h5file_masks.get_node(node_3, "/Z_{}/{}".format(z, sc))

    # Initialize the canvas
    canvas = numpy.zeros(shape=(width, height), dtype=numpy.float32)

    # Iterate all nodes: all trenches & rows
    for n in h5file_masks.walk_nodes(node_4, classname="CArray"):
        # Parse n's parent & its name to get keys
        row = n._v_parent._v_name
        trench = n.name.split("_")[1]
        mask = n.read()

        # Run a query on the data table & the mask, and return a modified float32 according to the function
        df_this_row_trench = df[
            (df["info_row_number"] == int(row))
            & (df["info_trench_number"] == int(trench))
        ]

        computed_image = numpy.zeros(shape=mask.shape, dtype=numpy.float32)
        for cell in df_this_row_trench.itertuples():
            mask_area = mask[
                cell.bounding_box_min_row : cell.bounding_box_max_row,
                cell.bounding_box_min_col : cell.bounding_box_max_col,
            ]

            computed_image[
                cell.bounding_box_min_row : cell.bounding_box_max_row,
                cell.bounding_box_min_col : cell.bounding_box_max_col,
            ] = (mask_area == cell.info_label) * eval(column_expression)

        # Write the result into the canvas
        key = "{}/region_{}".format(row, trench)
        r = regions[key]
        canvas[r[0] : r[2], r[1] : r[3]] = computed_image

    h5file_masks.close()

    return canvas.T


def computed_whole_labeled_canvas(
    file,
    fov,
    frame,
    z,
    width,
    height,
    sc,
    masks_file,
    column_expression,
    data_table_file,
):
    """
    FIXME Passing in width & height, even though they aren't needed, because it conforms
    to the generic shape shared by all the other functions. There must be a better way to do this!
    Maybe by using named arguments, or a dict of arguments?
    """
    h5file_masks = tables.open_file(masks_file, "r")
    # If no regions, then default to the zeroth region set number and region
    path = "/masks/{}/FOV_{}/Frame_{}/Z_{}/{}/0/region_0".format(
        file, fov, frame, z, sc
    )
    mask = h5file_masks.get_node(path).read()
    h5file_masks.close()

    # Run a query on the data table & the mask, and return a modified float32 according to the function
    df = query_for_masks(file, fov, frame, z, sc, data_table_file)

    # Speed up the compuation by only masking the bbox area.
    computed_image = numpy.zeros(shape=(width, height), dtype=numpy.float32)
    for cell in df.itertuples():
        mask_area = mask[
            cell.bounding_box_min_row : cell.bounding_box_max_row,
            cell.bounding_box_min_col : cell.bounding_box_max_col,
        ]

        computed_image[
            cell.bounding_box_min_row : cell.bounding_box_max_row,
            cell.bounding_box_min_col : cell.bounding_box_max_col,
        ] = (mask_area == cell.info_label) * eval(column_expression)

    return computed_image.T


### Image Layers


def add_image_layers(
    images_file,
    file_nodes,
    width,
    height,
    fields_of_view,
    frames,
    z_levels,
    viewer,
    channels,
    layer_params,
    corrections_dict,
):
    """
    Input an HDF5 file with images,
    image dimensions,
    fovs, frames, z dimensions,
    channel names,
    layer params.

    corrections_file: use to create a dict (input channel, output array for image to be divided by) [if no correction, set to 1.0]
    TODO: determine whether it's faster or simpler to set zero & one values for default, or to skip computation. For now, do it with zero and one.

    Adds all the image layers to the Napari viewer.
    """
    for c in channels:
        image_stack = dask_stack_loop(
            file_nodes,
            fields_of_view,
            frames,
            z_levels,
            width,
            height,
            load_img,
            numpy.float32,
            [numpy.float32, images_file, corrections_dict, c],
        )

        viewer.add_image(image_stack, **layer_params[c])


def add_computed_image_layers(
    masks_file,
    regions_file,
    file_nodes,
    fields_of_view,
    frames,
    z_levels,
    seg_channels,
    width,
    height,
    viewer,
    viewer_params,
    data_table_file,
):
    """
    Input an HDF5 file with images,
    image dimensions,
    fovs, frames, z dimensions,
    *segmentation* channel names,
    layer params.

    function: some sort of function which modifies the image

    data_table_file: HDF5 file containing a pytables table with measurements to use for computation

    Adds all the image layers to the Napari viewer.
    """
    if regions_file:
        for channel, categories in viewer_params.items():
            for category, items in categories.items():
                eval_string = items["eval_string"]
                params = items["viewer_params"]

                images_stack = dask_stack_loop(
                    file_nodes,
                    fields_of_view,
                    frames,
                    z_levels,
                    width,
                    height,
                    computed_masks_to_canvas,
                    numpy.float32,
                    [channel, regions_file, masks_file, eval_string, data_table_file],
                )

                viewer.add_image(images_stack, **params)

    # No regions file
    else:
        for channel, categories in viewer_params.items():
            for category, items in categories.items():
                eval_string = items["eval_string"]
                params = items["viewer_params"]

                images_stack = dask_stack_loop(
                    file_nodes,
                    fields_of_view,
                    frames,
                    z_levels,
                    width,
                    height,
                    computed_whole_labeled_canvas,
                    numpy.float32,
                    [channel, masks_file, eval_string, data_table_file],
                )

                viewer.add_image(images_stack, **params)


###


def get_largest_extents_hdf5(h5file, metadata_key):
    """
    TODO: for now, assume that the metadata types are integers,
    but in the future, could explicitly check, and set the min, max values accordingly.

    TODO: Also, support non list type objects.
    """
    from sys import maxsize

    smallest_value = maxsize
    largest_value = 0

    # only iterate the File_{} nodes
    metadata_node = h5file.get_node("/Metadata")()
    for node in h5file.iter_nodes(metadata_node):
        values = h5file.get_node(node, metadata_key).read()

        if values[0] < smallest_value:
            smallest_value = values[0]
        if values[-1] > largest_value:
            largest_value = values[-1]

    return (smallest_value, largest_value)


def metadata_attributes_equal(h5file, attribute):
    """
    Input an h5file with metadata, and an attribute of interest
    Returns None if not all nd2 files (in the h5file) have identical attributes of this type in their metadata
    Returns the attribute if they are identical
    """
    # only iterate the File_{} nodes
    metadata_node = h5file.get_node("/Metadata")()
    iter_nodes = h5file.iter_nodes(metadata_node)
    zeroth_attribute = h5file.get_node(next(iter_nodes), attribute).read()

    for node in iter_nodes:
        next_attribute = h5file.get_node(node, attribute).read()
        if next_attribute != zeroth_attribute:
            return None

    return zeroth_attribute


def metadata_array_equal(h5file, attribute):
    """
    Input an h5file with metadata, and an attribute of interest
    Returns None if not all nd2 files (in the h5file) have identical attributes of this type in their metadata
    Returns the attribute if they are identical
    TODO Update this desription. Is for comparing arrays!
    """
    # only iterate the File_{} nodes
    metadata_node = h5file.get_node("/Metadata")()
    iter_nodes = h5file.iter_nodes(metadata_node)
    zeroth_attribute = h5file.get_node(next(iter_nodes), attribute).read()

    for node in iter_nodes:
        next_attribute = h5file.get_node(node, attribute).read()
        if not numpy.array_equal(next_attribute, zeroth_attribute):
            return None

    return zeroth_attribute


def sub_bg_no_underflow(input_image, bg):
    """
    Subtract a background value from an image, but prevent underflow by stopping at zero.
    If input_image is a numpy array, then bg could be a fixed value,
    or another numpy array with the same dimensions as input_image.
    Dtype could be either floating point or integer.
    """
    return (input_image <= bg) * 0 + (input_image > bg) * (input_image - bg)


def load_img(
    file, fov, frame, z, width, height, dtype, images_file, corrections_dict, channel
):
    """
    Attempt to load an image, and if it cannot be found, return a zero'ed array instead.
    Apply camera bias and flat field corrections to the loaded image.
    """
    try:
        path = "/Images/{}/FOV_{}/Frame_{}/Z_{}/{}".format(file, fov, frame, z, channel)
        h5 = tables.open_file(images_file, "r")
        node = h5.get_node(path)
        img = node.read().astype(dtype)
        h5.close()

        if corrections_dict["camera_noise"] is not None:
            img = sub_bg_no_underflow(img, corrections_dict["camera_noise"])
        else:
            # A default value of 100.0
            img = sub_bg_no_underflow(img, 100.0)

        # Optional for each channel
        # If there is nothing, then do nothing!
        if corrections_dict["flatfield"] is not None:
            if corrections_dict["flatfield"][channel] is not None:
                img = numpy.divide(img, corrections_dict["flatfield"][channel])

        # Registration
        # Use try clause, because a None value will not be in the dict!
        if corrections_dict["registration"] is not None:
            registration_mat = corrections_dict["registration"][channel]
            if registration_mat is not None:
                print("Warping ", channel, "with matrix ", registration_mat)
                # FIXME it works on the entire 3D stack,
                # but we're only loading 1 Z level at a time.
                # --> How to apply the warping to a single slice?
                # Do we take the z-level, copy the frame into a zero'ed
                # arr at its slice level, apply the transform,
                # and then extract the slice again?
                img = imwarp(img, registration_mat, R_A=imref2d)

        return img.T

    except:
        return numpy.zeros((height, width), dtype=dtype)


def parse_corrections_settings(corrections_settings):
    """
    Input a dict with different possible corrections. Possible corrections include:
    - Camera noise (closed shutter)
    - Flatfielding (fluorescent dyes)
    - Fiducial registration (fluorescent beads)

    If specified, then load each of these into memory, and be able
    to pass it on through to the image layer loading.
    """
    # Initialize with default values, which might later be filled in.
    corrections_dict = {"camera_noise": None, "flatfield": None, "registration": None}

    # Was a corrections file specified?
    # If not, then just pass back some Nones
    try:
        corr_h5 = tables.open_file(corrections_settings["corrections_file"])

        # Camera noise
        try:
            camera_noise = corr_h5.get_node("/CAMERA/camera_noise").read()
            corrections_dict["camera_noise"] = camera_noise
        except:
            corrections_dict["camera_noise"] = None

        # Flat-fielding matrices
        try:
            ff_dict = {}
            for ch_img, ch_corr in corrections_settings[
                "flatfield_channel_mapping"
            ].items():
                # Optional for each channel
                if ch_corr:
                    corr_img = corr_h5.get_node("/FLATFIELD/{}".format(ch_corr)).read()
                    ff_dict[ch_img] = corr_img
                else:
                    ff_dict[ch_img] = 1.0

            corrections_dict["flatfield"] = ff_dict

        # Default value of 1.0 (no change)
        except:
            corrections_dict["flatfield"] = None

        # Registration calibration (fiducial beads, e.g. TetraSpeck)
        try:
            # Input reference_moving_mode_inworld, output array
            registration_dict = {}

            # Open the HDF5 file with registration matrices
            registration_h5 = tables.open_file(
                corrections_settings["registration"]["matrix_file"], "r"
            )

            for ch, vals in corrections_settings["registration"]["mapping"].items():
                # None is used for the reference channel in the settings file
                if vals is not None:
                    path = "/{}/{}/{}/inworld_{}".format(
                        vals["reference"], vals["moving"], vals["mode"], vals["inworld"]
                    )
                    registration_dict[ch] = registration_h5.get_node(path).read()
                # None is used for the reference channel in the registrations dict
                else:
                    registration_dict[ch] = None

            registration_h5.close()
            corrections_dict["registration"] = registration_dict

        except:
            corrections_dict["registration"] = None

        # Done!
        corr_h5.close()

    except:
        print("Warning: could not open corrections file!")

    return corrections_dict


### Main function


def main_hdf5_browser_function(
    images_file,
    masks_file,
    regions_file,
    settings_file,
    viewer_params_file,
    data_table_file,
):
    """
    Use Napari to browse an HDF5 file with microscope images arranged in a hierarchy:
    ND2 File -> Field of View -> Time Frame -> Z_Level -> Channel
    Images are loaded on demand through use of a "lazy" dask array.

    If specified, can also browse labeled regions, corresponding to either trenches or segmented cells.
    Use a YAML file to specify the Napari settings.
    """
    # Initialilize the viewer
    viewer = napari.Viewer()

    # Global settings (YAML)
    settings = read_params_file(settings_file)

    # Load the image layer params for napari
    # TODO params for masks layers? regions layers?
    layer_params = settings["napari_settings"]

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

    # Load correction settings
    # For camera noise, flatfielding, registration
    try:
        corrections_settings = settings["corrections"]
    except:
        corrections_settings = None

    corrections_dict = parse_corrections_settings(corrections_settings)

    # Image layers
    add_image_layers(
        images_file,
        file_nodes,
        width,
        height,
        extents["fields_of_view"],
        extents["frames"],
        extents["z_levels"],
        viewer,
        channels,
        layer_params,
        corrections_dict,
    )

    # Regions (trench) layer
    # TODO call this multiple times, if there are multiple trench rows?
    # (With each trench row getting its own layer)
    # Requires querying the table & determining the number of unique trench rows.
    if regions_file:
        add_regions_layer(
            regions_file,
            file_nodes,
            extents["fields_of_view"],
            extents["frames"],
            extents["z_levels"],
            viewer,
            width,
            height,
        )

    # Masks layers
    if masks_file:
        masks_file_path = os.path.join(masks_file, "MASKS/masks.h5")
        h5file_masks = tables.open_file(masks_file_path, "r")

        # Load the segmentation channels from the masks h5file
        seg_params_node = h5file_masks.get_node("/Parameters", "seg_params.yaml")
        seg_params_str = seg_params_node.read().tostring().decode("utf-8")
        seg_params_dict = read_params_string(seg_params_str)
        channels_list = seg_params_dict["segmentation"]
        seg_channels = [k for k in channels_list.keys()]

        h5file_masks.close()

        add_masks_layer(
            masks_file_path,
            regions_file,
            file_nodes,
            extents["fields_of_view"],
            extents["frames"],
            extents["z_levels"],
            seg_channels,
            width,
            height,
            viewer,
        )

    # Would it make more sense to pre-compute the columns? Basically
    # make it a column display, and just pass in the name of the column.
    # Input a set of images, but also run a calculation on them before displaying.
    if viewer_params_file and data_table_file and masks_file:
        viewer_params = read_params_file(viewer_params_file)
        masks_file_path = os.path.join(masks_file, "MASKS/masks.h5")

        add_computed_image_layers(
            masks_file_path,
            regions_file,
            file_nodes,
            extents["fields_of_view"],
            extents["frames"],
            extents["z_levels"],
            seg_channels,
            width,
            height,
            viewer,
            viewer_params,
            data_table_file,
        )

    # Launch the viewer
    napari.run()
