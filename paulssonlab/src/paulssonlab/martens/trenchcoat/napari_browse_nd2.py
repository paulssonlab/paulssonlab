import os

import dask.array
import napari
import nd2reader
import numpy
from convert import make_nd2_list
from dask import delayed
from metadata import get_largest_extents, metadata_attributes_equal
from params import read_params_file

"""
Use Napari to browse a directory of ND2 files.

TODO: what if each "file" within the ND2 directory had different channels?
For example, this comes up when flatfielding using dyes. Each dye has only a subset of the channels.
There should be a way to insert zeros as needed in the empty positions.
"""


def main_nd2_browser_function(in_dir, settings_file):
    """
    Open a directory of ND2 files in Napari, using the parameters specified in the YAML file.
    """
    trenchcoat_settings = read_params_file(settings_file)

    # 1. Load the params, and use the params to determine the list of channel names
    layer_params = trenchcoat_settings["napari_settings"]
    channels = layer_params.keys()

    # 2. Open up all the nd2 files in a dict., and determine the max. ranges of FOVs, frames, z_levels across all files.
    nd2_files = make_nd2_list(in_dir)

    nd2_dict = {}
    for f in nd2_files:
        nd2_dict[f] = nd2reader.Nd2(f)

    attributes = ["fields_of_view", "frames", "z_levels"]
    extents = {}
    for a in attributes:
        (smallest, largest) = get_largest_extents(nd2_dict, a)
        if smallest == largest:
            extents[a] = [smallest]
        else:
            extents[a] = [i for i in range(smallest, largest + 1)]

    height = metadata_attributes_equal(nd2_dict, "height")
    width = metadata_attributes_equal(nd2_dict, "width")

    # 3. Init napari
    viewer = napari.Viewer()

    # 4. Load the images into a dask array using delayed (lazy) loading
    # define a lazy function once, for later use
    delayed_zeros = delayed(numpy.zeros)

    # NOTE: as written this seems awfully confusing, with so many loops and dicts and lists,
    # but it makes sense. The idea is that each set of images for a channel, in a given layer,
    # gets loaded up into a list, which supports appending, then is converted to a dask array
    # with lazy loading.
    #
    # In order to keep track of each channel's dask array, it has to be stored into
    # a dict. Then, at each level of stacking, the same process is repeated.
    #
    # Finally, we have a large dask array for each channel, spanning all the FOVs, Frames, and Z_levels.
    #
    # The fundamental reason for doing it this way, instead of looping the channels on the outside,
    # is that the nd2 file format stores the images interspersed with one another. This means that
    # loading an image stack with many channels is progressively more wasteful, as a bunch of data
    # are loaded from disk and then discarded!
    #
    # By loading the data at once, this expensive task is done just once, and instead the cheap looping is repeated.
    #
    # Also note: https://github.com/dask/dask/issues/2000 Cannot modify dask arrays in-place...
    # This means we can't use a much more sane way of iterating + filling in with delayed data. (Yes, I did try it.)

    # Files
    file_array = {}
    for c in channels:
        file_array[c] = []

    # TODO sorted iteration order?
    for _, reader in nd2_dict.items():
        delayed_get_image_stack = delayed(reader.get_image_stack)
        # FOVs
        fov_array = {}
        for c in channels:
            fov_array[c] = []

        for fov in extents["fields_of_view"]:
            # Frames
            frame_array = {}
            for c in channels:
                frame_array[c] = []

            for frame in extents["frames"]:
                # Z levels
                z_array = {}
                for c in channels:
                    z_array[c] = []

                for z in extents["z_levels"]:
                    # 1. load stack of images from all channels, given this fov, frame, z_level
                    # NOTE: this uses the modified function (not present in original nd2reader library)

                    # a. if images were found: iterate channels & add images to the dask arrays, lazily
                    try:
                        stack = delayed_get_image_stack(
                            field_of_view=fov, frame_number=frame, z_level=z
                        )

                    # b. if images cannot be found, then for each channel add zeros array with same dimension (lazily)
                    except:
                        # FIXME width, height -- or height, width??
                        stack = delayed_zeros(
                            shape=(width, height, len(channels)),
                            dtype=numpy.uint16,
                            order="F",
                        )

                    # We enumerate reader.channels so that the index i matches the channel c.
                    for i, c in enumerate(reader.channels):
                        z_array[c].append(
                            dask.array.from_delayed(
                                stack[..., i],
                                shape=(width, height),
                                dtype=numpy.uint16,
                            )
                        )

                for c in channels:
                    z_array_to_stack = dask.array.stack(z_array[c])
                    frame_array[c].append(z_array_to_stack)

            for c in channels:
                frame_array_to_stack = dask.array.stack(frame_array[c])
                fov_array[c].append(frame_array_to_stack)

        for c in channels:
            fov_array_to_stack = dask.array.stack(fov_array[c])
            file_array[c].append(fov_array_to_stack)

    for c in channels:
        megastack = dask.array.stack(file_array[c])
        viewer.add_image(megastack, **layer_params[c])

    napari.run()
