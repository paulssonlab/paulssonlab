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

### Regions (e.g. Trenches) Layer


def add_regions_layer(regions_table, viewer):
    """
    There are 2 kinds of regions:
    1. Trench rows
    2. Trenches

    A dataset could have either one, the other, or both.
    In each case, the individual elements (rows, trenches) must be pasted into a canvas which spans the original image.

    Create one label layer for trench rows, and another label layer for trenches.
    TODO The trenches from different rows will be placed within the same, or different, labels layers???
    Doesn't it make the most sense for each to have its own layer? why not?

    TODO pass in information for the labels layer params? name etc.

    FIXME what if the regions are shared across Z, and therefore the dimensionality of this labels layer
    is not the same as the dims of the images or the masks? Will it "just work?"
    """

    # TEMP: just do 1 row of trenches, and not the row itself.
    trenches_stack = make_regions_lazy_dask_stack(regions_table)
    viewer.add_labels(trenches_stack)  # , **layer_params) # TODO layer params?

    ## 1. Trench rows
    # trench_rows_stack = make_regions_lazy_dask_stack(h5file, file_nodes, fields_of_view, frames, "row_coords")
    # viewer.add_labels(trench_rows_stack) #, **layer_params) # TODO layer params?

    # row_coords_arr = h5file.get_node().read()
    # num_trench_rows = row_coords_arr.shape[0]

    ## FIXME what if some FOVs have 1 trench row, and others have 2 trench rows? How to handle different numbers
    ## in the Dask array? Need to determine the largest possible number of rows, and then fill in
    ## 2. Trenches within each row
    # for n in num_trench_rows:
    # trenches_stack = make_regions_lazy_dask_stack(h5file, file_nodes, fields_of_view, frames, "row_coords/row_{}_trench_coords".format(n))
    # viewer.add_labels(trenches_stack) #, **layer_params) # TODO layer params?


def make_regions_lazy_dask_stack(regions_table):
    """
    Input the location of the region (trench row or trenches) coordinates within the HDF5 file,
    output a dask array stacked for all the dimensions,
    whose final dimension is the canvas with the labeled regions,
    representing either trench rows, or trenches within a row.
    """
    # Compile a canvas with all the trench regions
    lazy_regions_arr_to_canvas = delayed(regions_arr_to_canvas)

    whole_table = regions_table.read()
    whole_df = pandas.DataFrame(whole_table)
    file_nodes = whole_df["info_file"].unique()
    fields_of_view = whole_df["info_fov"].unique()
    frames = whole_df["info_frame"].unique()

    file_array = []
    for file in file_nodes:
        # FOV nodes
        fov_array = []
        for fov in fields_of_view:
            # Frame nodes
            frame_array = []
            for frame in frames:
                # Regions are always shared across Z, so it's OK to read them in at this level.
                # path = "/{}/FOV_{}/Frame_{}/{}".format(file, fov, frame, node_path)
                # FIXME pass in image width & height
                arr = dask.array.from_delayed(
                    lazy_regions_arr_to_canvas(
                        regions_table, file, fov, frame, 2048, 2048
                    ),
                    shape=(2048, 2048),
                    dtype=numpy.uint16,
                )
                frame_array.append(arr)

            frame_array_to_stack = dask.array.stack(frame_array)
            fov_array.append(frame_array_to_stack)

        fov_array_to_stack = dask.array.stack(fov_array)
        file_array.append(fov_array_to_stack)

    megastack = dask.array.stack(file_array)

    return megastack


def regions_arr_to_canvas(table, file, fov, frame, width, height):
    """
    ## Input a numpy array with many rows, each containing the rectangular coordinates of regions.

    Pass in info from which to read corresponding trench coords from the global trench coords table.

    Return a numpy array with numbered pixel positions for each region, and zeros for background.
    NOTE: what happens if the regions overlap?
    """
    # regions_coords = h5file_regions.get_node(path).read()
    # TODO: trench row number?
    condition = "(info_file == file) & (info_fov == fov) & (info_frame == frame)"
    search = table.read_where(condition)
    df = pandas.DataFrame(search)

    regions = []
    for r in df.itertuples():
        regions.append([r.min_row, r.min_col, r.max_row, r.max_col])

    regions = numpy.array(regions)

    # FIXME or height, width?
    canvas = numpy.zeros(shape=(width, height), dtype=numpy.uint16)

    for i, r in enumerate(regions):
        canvas[r[0] : r[2], r[1] : r[3]] = i + 1

    return numpy.rot90(canvas, k=-1)


### Masks (e.g. cells) Layer


def add_masks_layer(
    masks_file, regions_file, fields_of_view, frames, z_levels, height, width, viewer
):
    """
    Because each region gets its own masks, they need to be merged!
    This means that it is necessary to specify which regions were used to do so during the segmentation step.

    TODO figure out how to require regions? (or, the same region will just be over-written, and the user will figure it out?)

    As with the trenches, these regions the masks within these regions need to be combined into a single, large image.

    TODO pass in information for the labels layer params? name etc.
    """
    h5file_masks = tables.open_file(masks_file, "r")

    # Load the segmentation channels from the masks h5file
    seg_params_node = h5file_masks.get_node("/Parameters", "seg_params.yaml")
    seg_params_str = seg_params_node.read().tostring().decode("utf-8")
    seg_params_dict = read_params_string(seg_params_str)
    seg_channels = [k for k in seg_params_dict.keys()]

    file_nodes = [x._v_name for x in h5file_masks.list_nodes("/masks")]

    # If there were regions, then the masks need to be agglomerated into a larger image
    if regions_file:
        # Compile a canvas with all the masks from all the segmented areas
        lazy_masks_to_canvas = delayed(masks_to_canvas)

        h5file_regions = tables.open_file(regions_file, "r")
        table = h5file_regions.root.trench_coords

        # FIXME how to determine the number of region sets? should each region set get its own labels layer?
        # that would mean sc x region_sets total number of layers.
        #
        # NOTE that if the regions within a set overlap, then they will overwrite each other
        # when it comes time to merging them together in the same canvas.
        for sc in seg_channels:
            # File nodes
            # FIXME what order to iterate? do they need to be sorted?
            # Needs to be the same order as the images were loaded.
            file_array = []
            for file in file_nodes:
                # Go through series of de-references b/c pytables doesn't allow iterating
                # over nodes if some of the nodes are symlinks.
                path = "/masks/{}".format(file)
                node_1 = h5file_masks.get_node(path)()

                # FOV nodes
                fov_array = []
                for fov in fields_of_view:
                    node_2 = h5file_masks.get_node(node_1, "/FOV_{}".format(fov))()

                    # NOTE: for now, trenches are shared across frames
                    # TODO: what if they aren't?
                    condition = "(info_file == file) & (info_fov == fov)"
                    search = table.read_where(condition)
                    df = pandas.DataFrame(search)

                    # Frame nodes
                    frame_array = []
                    for frame in frames:
                        node_3 = h5file_masks.get_node(
                            node_2, "/Frame_{}".format(frame)
                        )()

                        # Regions are always shared across Z, so it's OK to read them in now.
                        # FIXME read them in from the table.
                        # path = "/{}/FOV_{}/Frame_{}"
                        # regions_coords = h5file_regions.get_node(path).read()

                        regions = {}
                        for r in df.itertuples():
                            key = "{}/region_{}".format(
                                r.info_row_number, r.info_trench_number
                            )
                            regions[key] = [r.min_row, r.min_col, r.max_row, r.max_col]

                        # Z nodes
                        z_array = []
                        for z in z_levels:
                            # TODO now, need to loop all the available regions,
                            # and then copy over the corresponding data from the masks into
                            # a single, "zeroed out" canvas
                            # FIXME issue iterating nodes when accessing via symlink!
                            # this is bad, because that would require opening separate files
                            # depending on the data, but how do we lazy load from separate files?
                            node_4 = h5file_masks.get_node(
                                node_3, "/Z_{}/{}".format(z, sc)
                            )

                            # FIXME how do we know that the nodes iteration is the same order as the numbering of the trenches?
                            # could be non-numeric sorting)
                            # FIXME shape
                            arr = dask.array.from_delayed(
                                lazy_masks_to_canvas(
                                    regions, h5file_masks, node_4, width, height
                                ),
                                shape=(2048, 2048),
                                dtype=numpy.uint16,
                            )
                            z_array.append(arr)

                        z_array_to_stack = dask.array.stack(z_array)
                        frame_array.append(z_array_to_stack)

                    frame_array_to_stack = dask.array.stack(frame_array)
                    fov_array.append(frame_array_to_stack)

                fov_array_to_stack = dask.array.stack(fov_array)
                file_array.append(fov_array_to_stack)

            # TODO what params to pass in? use a yaml file?
            megastack = dask.array.stack(file_array)
            viewer.add_labels(megastack)

    ### For reference:
    ### h5file_masks.create_carray("/{}/{}/{}".format(z_level, sc, region_set_number),

    # No regions file
    else:
        for sc in seg_channels:
            # File nodes
            # FIXME what order to iterate? do they need to be sorted?
            # Needs to be the same order as the images were loaded.
            file_array = []
            for file in file_nodes:
                # FOV nodes
                fov_array = []
                for fov in fields_of_view:
                    # Frame nodes
                    frame_array = []
                    for frame in frames:
                        # Z nodes
                        z_array = []
                        for z in z_levels:
                            # If no regions, then default to the zeroth region set number and region
                            path = (
                                "/masks/{}/FOV_{}/Frame_{}/Z_{}/{}/0/region_0".format(
                                    file, fov, frame, z, sc
                                )
                            )
                            node = h5file_masks.get_node(path)
                            # Rotate by -90 deg to account for F -> C ordering of arrays when visualizing in Napari
                            arr = dask.array.from_delayed(
                                delayed(numpy.rot90(node.read(), k=-1)),
                                # If there are no regions, then it must fill the entire
                                # dimensions of the image.
                                shape=(height, width),  # FIXME or width,height?
                                dtype=numpy.uint16,
                            )

                            z_array.append(arr)

                        z_array_to_stack = dask.array.stack(z_array)
                        frame_array.append(z_array_to_stack)

                    frame_array_to_stack = dask.array.stack(frame_array)
                    fov_array.append(frame_array_to_stack)

                fov_array_to_stack = dask.array.stack(fov_array)
                file_array.append(fov_array_to_stack)

            # TODO what params to pass in? use a yaml file?
            megastack = dask.array.stack(file_array)
            viewer.add_labels(megastack)


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


def masks_to_canvas(regions, h5file_masks, node, width, height):
    """
    ###Input a numpy array with many rows, each containing the rectangular coordinates of regions.
    Input a dict with keys structured as "rownum/region_regionnum", each containing the rectangular coordinates of regions.

    Input an iterator over an H5 file node, below which hang individual segmentation masks.
    Return a numpy array with labeled pixel positions across all the masks, and zeros for background.
    NOTE: what happens if the regions overlap?
    """
    # FIXME or height, width?
    canvas = numpy.zeros(shape=(width, height), dtype=numpy.uint16)

    for n in h5file_masks.walk_nodes(node, classname="CArray"):
        # Parse n's parent & its name to get keys
        row = n._v_parent._v_name
        trench = n.name.split("_")[1]
        key = "{}/region_{}".format(row, trench)
        r = regions[key]
        canvas[r[0] : r[2], r[1] : r[3]] = n.read()

    return numpy.rot90(canvas, k=-1)


### Image Layers


def add_image_layers(
    h5file,
    height,
    width,
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
    # Open & parse the corrections file
    corrections = None
    if corrections_file:
        corrections = load_values(corrections_file)

    # Define a function for lazily loading a single image from the h5 file
    lazy_load_img = delayed(load_img)

    # Iterate the H5 file images & lazily load them into a dask array
    file_nodes = [x._v_name for x in h5file.list_nodes("/Images")]

    for c in channels:
        c = c.decode("utf-8")

        # If they weren't defined, then define them here
        # TODO allow specifying the name of the channel?
        try:
            bias = corrections["DARK"]
        # No camera bias correction, so subtract zero.
        except:
            bias = 0.0

        try:
            ffc = corrections[c]
        # No flat-field inhomogeneity correction, so divide by 1.0.
        except:
            ffc = 1.0

        # File nodes
        # FIXME what order to iterate? do they need to be sorted?
        file_array = []
        for file in file_nodes:
            # FOV nodes
            fov_array = []
            for fov in fields_of_view:
                # Frame nodes
                frame_array = []
                for frame in frames:
                    # Z nodes
                    z_array = []
                    for z in z_levels:
                        path = "/Images/{}/FOV_{}/Frame_{}/Z_{}".format(
                            file, fov, frame, z
                        )
                        node = h5file.get_node(path)
                        # FIXME problem here: what if camera_biases is None? Then we can't index it!
                        arr = dask.array.from_delayed(
                            lazy_load_img(
                                node, c, height, width, numpy.float32, bias, ffc
                            ),
                            shape=(height, width),  # FIXME or width,height?
                            dtype=numpy.float32,
                        )
                        z_array.append(arr)

                    z_array_to_stack = dask.array.stack(z_array)
                    frame_array.append(z_array_to_stack)

                frame_array_to_stack = dask.array.stack(frame_array)
                fov_array.append(frame_array_to_stack)

            fov_array_to_stack = dask.array.stack(fov_array)
            file_array.append(fov_array_to_stack)

        megastack = dask.array.stack(file_array)
        viewer.add_image(megastack, **layer_params[c])


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


def load_img(node, channel, height, width, dtype, camera_bias, flatfield_correction):
    """
    Attempt to load an image, and if it cannot be found, return a zero'ed array instead.
    Apply camera bias and flat field corrections.
    We use Fortran indexing, Napari usese C indexing, so rotate the image by -90 degrees.
    """
    try:
        # return (node._f_get_child(channel).read() - camera_bias) / flatfield_correction
        return numpy.rot90(
            (node._f_get_child(channel).read() - camera_bias) / flatfield_correction,
            k=-1,
        )

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

        # Image layers
        add_image_layers(
            h5file,
            height,
            width,
            extents["fields_of_view"],
            extents["frames"],
            extents["z_levels"],
            viewer,
            channels,
            layer_params,
            corrections_file,
        )

        # Masks layers
        # FIXME issue with opening/closing regions h5 file & lazy loading of the data...
        if masks_file:
            add_masks_layer(
                masks_file,
                regions_file,
                extents["fields_of_view"],
                extents["frames"],
                extents["z_levels"],
                height,
                width,
                viewer,
            )

        # Regions (trench) layer
        if regions_file:
            h5file_regions = tables.open_file(regions_file, "r")
            regions_table = h5file_regions.root.trench_coords

            add_regions_layer(
                regions_table,
                # extents["fields_of_view"],
                # extents["frames"],
                viewer,
            )

    # Done!
    h5file.close()

    if regions_file:
        h5file_regions.close()
