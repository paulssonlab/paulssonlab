#!/usr/bin/env python3

import tables
import numpy

from metadata import get_metadata
from params import read_params_file


def main_corrections_function(in_file, out_file, dark_channel, bg_file):
    """
    Input path to an H5 file which has images used to correct two aspects of microscope images:
    1. Camera bias: per-pixel intensities which are non-zero
    2. Uneven illumination: per-pixel intensities which are brighter or dimmer depending on how close they are to the edges.

    Input path for writing the matrices to a new H5 file.

    Also input either the name of a dark channel, OR a dict containing single background values per channel.

    By convention, any images in a "DARK" channel will be assumed to correspond to images taken without any illumination.
    These are used to fix camera bias. All other images will have the DARK image subtracted from them.
    In the absence of a DARK channel, another option is to subtract a fixed value. This can still be helpful in case
    the camera artifically adds some value (e.g. 100) to every pixel position, to prevent negative intensities
    where the bias tends towards less than the mean bias.

    It is assumed that pixel values never are less than 0.0. All images are converted to float32.

    The general procedure is to:
    - compute the mean of all images of a given type,
    - subtract the bias
    - normalize by the mean pixel intensity.

    Normalization by the mean intensity is chosen over the max intensity as a way to accomodate hot pixels. In this way, these normalization
    matrices will tend to be more reproducible, and to tend towards values close to 1.0.

    The alternative convention would always map values between 0.0 and 1.0, which has the advantage of preserving intensity ranges.

    Should the preference for mean vs max change, this can easily be done by modifying "arr /= arr.mean()" to be "arr /= arr.max()".
    NOTE: could make this an option.

    NOTE: could also add options such as running local window averaging, to smooth out variability. Would also help with hot pixels.
    """

    if dark_channel and bg_file:
        # FIXME Error!
        pass
    elif bg_file:
        background_values = read_params_file(bg_file)

    h5file = tables.open_file(in_file, mode="r")

    file_names = [i._v_name for i in h5file.list_nodes("/Images")]

    node = h5file.get_node("/Metadata/{}".format(file_names[0]))()
    channels = get_metadata(node)["channels"]
    width = get_metadata(node)["width"]
    height = get_metadata(node)["height"]

    total = 0
    for n in h5file.iter_nodes(h5file.root.Images):
        metadata_node = h5file.get_node("/Metadata/{}".format(n._v_name))()
        metadata = get_metadata(metadata_node)
        total += (
            len(metadata["fields_of_view"])
            * len(metadata["frames"])
            * len(metadata["z_levels"])
        )

    ### Calculate the means
    # Input a channel name, output a numpy array which is the average of many images
    means = {}

    for c in channels:
        c = c.decode("utf-8")
        # FIXME Or is it width, height?
        stack = numpy.empty(
            shape=(height, width, total), dtype=numpy.float32, order="F"
        )

        i = 0
        for n in h5file.iter_nodes(h5file.root.Images):

            metadata_node = h5file.get_node("/Metadata/{}".format(n._v_name))()
            metadata = get_metadata(metadata_node)

            for fov in metadata["fields_of_view"]:
                for frame in metadata["frames"]:
                    for z_level in metadata["z_levels"]:
                        path = "/Images/{}/FOV_{}/Frame_{}/Z_{}/{}".format(
                            n._v_name, fov, frame, z_level, c
                        )
                        img_node = h5file.get_node(path)
                        img = img_node.read()
                        stack[..., i] = img
                        i += 1

        # Compute the mean
        # TODO / FIXME: would it be more memory efficient to sum then manually, and then divide by the number?
        # Would this come at the expense of precision? (not sure how numpy handles the mean calculation)
        # Don't have to worry about overflow, because we are only working with float32.
        means[c] = stack.mean(axis=2)

    h5file.close()

    ### Normalize & write to file
    corrections_h5 = tables.open_file(out_file, mode="w")

    for c, arr in means.items():
        # For the DARK channel, just write it without further normalization (preserve the absolute intensities)
        if c == dark_channel:
            corrections_h5.create_carray(
                "/",
                c,
                obj=arr,
                chunkshape=(128, 128),
                # Complevel 9 vs 1 is not worth it. Very minimal gains in compression amount.
                filters=tables.Filters(complevel=1, complib="zlib"),
            )

        # Skip these steps for DARK channel (regardless of whether not it exists)
        else:
            # Check whether DARK channel actually exists
            if dark_channel:
                # If it does, subtract it from the array
                arr -= means[dark_channel]
            else:
                # No dark channel, so option to subtract a constant value instead
                # arr -= background_values[c]
                bg_value = background_values[c]
                # Don't allow the values to wrap around
                # NOTE bg_value could also be an array?
                arr = (arr <= bg_value) * 0.0 + (arr > bg_value) * (arr - bg_value)

            # Divide by the mean value to normalize all values
            # Dividing by the mean, as opposed to the max, reduces the effect of outlier "hot" pixels
            arr /= arr.mean()

            corrections_h5.create_carray(
                "/",
                c,
                obj=arr,
                chunkshape=(128, 128),
                # Complevel 9 vs 1 is not worth it. Very minimal gains in compression amount.
                filters=tables.Filters(complevel=1, complib="zlib"),
            )

    corrections_h5.close()
