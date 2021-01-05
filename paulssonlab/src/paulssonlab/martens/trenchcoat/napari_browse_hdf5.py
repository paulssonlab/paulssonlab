import tables
import numpy
import napari
from dask import delayed
import dask.array
import pandas
from params import read_params_string, read_params_file

"""
NOTE for now, assumes image dtype is float32

TODO clean up the metadata functions, and their names, then stick them in the separate metadata Python file
Might also need to rework the similar functions used for nd2 browsing?
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
    query = "(info_file == file) & (info_fov == fov) & (info_frame == frame) & (info_z_level == z_level)"
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

    # FIXME or height, width?
    canvas = numpy.zeros(shape=(width, height), dtype=numpy.uint16)

    for r, tr in zip(regions, trench_numbers):
        # +1 because zero is background
        canvas[r[0] : r[2], r[1] : r[3]] = tr + 1

    return canvas.T


def regions_arr_to_canvas(file, fov, frame, z_level, width, height, regions_file):
    """
    Input a numpy array with many rows, each containing the rectangular coordinates of regions.

    NOTE now, the input is a table containing relevant data,
    and then the table is queried to create the regions arr.

    Pass in info from which to read corresponding trench coords from the global trench coords table.

    Return a numpy array with numbered pixel positions for each region, and zeros for background.
    NOTE: what happens if the regions overlap?
    """
    # TODO: trench row number?
    h5file_regions = tables.open_file(regions_file, "r")
    regions_table = h5file_regions.get_node("/", "row_coords")

    condition = "(info_file == file) & (info_fov == fov) & (info_frame == frame) & (info_z_level == z_level)"
    search = regions_table.read_where(condition)
    df = pandas.DataFrame(search)
    h5file_regions.close()

    regions = []
    for r in df.itertuples():
        regions.append([r.min_row, r.min_col, r.max_row, r.max_col])

    # FIXME or height, width?
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

    # FIXME how to determine the number of region sets? should each region set get its own labels layer?
    # that would mean sc x region_sets total number of layers.
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
    # FIXME need to re-work the lazy part! Follow the lead from the others.
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


def load_values(in_file):
    """
    Load values from an HDF5 file into a dict. Useful for flatfield corrections & camera biases.
    TODO find a better name for this subroutine
    """
    if in_file:
        dict = {}

        h5 = tables.open_file(in_file, "r")
        for node in h5.iter_nodes("/"):
            key = node._v_name
            dict[key] = node.read()
        h5.close()

    else:
        dict = None

    return dict


def whole_labeled_canvas(file, fov, frame, z, width, height, sc, masks_file):
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
    img = h5file_masks.get_node(path).read()
    h5file_masks.close()

    # Rotate by -90 deg to account for F -> C ordering of arrays when visualizing in Napari
    return numpy.rot90(img, k=-1)


def masks_to_canvas(file, fov, frame, z, width, height, sc, regions_file, masks_file):
    """
    TODO update documentation how this function works
    Input a numpy array with many rows, each containing the rectangular coordinates of regions.

    Input a dict with keys structured as "rownum/region_regionnum", each containing the rectangular coordinates of regions.

    Return a numpy array with labeled pixel positions across all the masks, and zeros for background.
    NOTE: what happens if the regions overlap?
    """
    h5file_regions = tables.open_file(regions_file, "r")
    table_regions = h5file_regions.get_node("/trench_coords")
    condition = "(info_file == file) & (info_fov == fov) & (info_frame == frame) & (info_z_level == z)"
    search_regions = table_regions.read_where(condition)
    h5file_regions.close()
    df_regions = pandas.DataFrame(search_regions)

    # Create a dict: input row number & trench number,
    # output bbox coords of the corresponding trench.
    regions = {}
    for r in df_regions.itertuples():
        key = "{}/region_{}".format(r.info_row_number, r.info_trench_number)
        regions[key] = (r.min_row, r.min_col, r.max_row, r.max_col)

    h5file_masks = tables.open_file(masks_file, "r")

    # Go through series of de-references b/c pytables doesn't allow iterating
    # over nodes if some of the nodes are symlinks.
    path = "/masks/{}".format(file)
    node_1 = h5file_masks.get_node(path)()
    node_2 = h5file_masks.get_node(node_1, "/FOV_{}".format(fov))()
    node_3 = h5file_masks.get_node(node_2, "/Frame_{}".format(frame))()
    node_4 = h5file_masks.get_node(node_3, "/Z_{}/{}".format(z, sc))

    # FIXME or height, width?
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
    corrections_file,
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
    # Define a function for lazily loading a single image from the h5 file
    lazy_load_img = delayed(load_img)

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
            [numpy.float32, images_file, corrections_file, c],
        )

        viewer.add_image(image_stack, **layer_params[c])


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


def load_img(
    file, fov, frame, z, width, height, dtype, images_file, corrections_file, channel
):
    """
    Attempt to load an image, and if it cannot be found, return a zero'ed array instead.
    Apply camera bias and flat field corrections.
    We use Fortran indexing, Napari usese C indexing, so rotate the image by -90 degrees.
    """
    try:
        # Open & parse the corrections file
        corrections = None
        if corrections_file:
            corrections = load_values(corrections_file)

        # If they weren't defined, then define them here
        # TODO allow specifying the name of the channel?
        try:
            camera_bias = corrections["DARK"]
        # No camera bias correction, so subtract zero.
        except:
            camera_bias = 0.0

        try:
            flatfield_correction = corrections[channel]
        # No flat-field inhomogeneity correction, so divide by 1.0.
        except:
            flatfield_correction = 1.0

        path = "/Images/{}/FOV_{}/Frame_{}/Z_{}".format(file, fov, frame, z)
        h5 = tables.open_file(images_file, "r")
        node = h5.get_node(path)
        img = (node._f_get_child(channel).read() - camera_bias) / flatfield_correction
        h5.close()
        return img.T

    except:
        # FIXME width,height?
        return dask.array.from_delayed(numpy.zeros(shape=(height, width), dtype=dtype))


### Main function


def main_hdf5_browser_function(
    images_file, masks_file, regions_file, corrections_file, napari_settings_file
):
    """
    Use Napari to browse an HDF5 file with microscope images arranged in a hierarchy:
    ND2 File -> Field of View -> Time Frame -> Z_Level -> Channel
    Images are loaded on demand through use of a "lazy" dask array.

    If specified, can also browse labeled regions, corresponding to either trenches or segmented cells.
    Use a YAML file to specify the Napari settings.
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
            corrections_file,
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
            h5file_masks = tables.open_file(masks_file, "r")

            # Load the segmentation channels from the masks h5file
            seg_params_node = h5file_masks.get_node("/Parameters", "seg_params.yaml")
            seg_params_str = seg_params_node.read().tostring().decode("utf-8")
            seg_params_dict = read_params_string(seg_params_str)
            seg_channels = [k for k in seg_params_dict.keys()]

            h5file_masks.close()

            add_masks_layer(
                masks_file,
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
