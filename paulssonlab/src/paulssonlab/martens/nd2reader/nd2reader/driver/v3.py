# -*- coding: utf-8 -*-

import array
import numpy as np
import struct
from nd2reader.model.image import Image
from nd2reader.common.v3 import read_chunk
from nd2reader.exc import NoImageError


class V3Driver(object):
    """
    Accesses images from ND2 files made with NIS Elements 4.x. Confusingly, files of this type have a version number of 3.0+.

    """

    def __init__(self, metadata, label_map, file_handle):
        """
        :param metadata:    a Metadata object
        :param label_map:    a raw dictionary of pointers to image locations
        :param file_handle:    an open file handle to the ND2

        """
        self._metadata = metadata
        self._label_map = label_map
        self._file_handle = file_handle

    def calculate_image_properties(self, index):
        field_of_view = self._calculate_field_of_view(index)
        channel = self._calculate_channel(index)
        z_level = self._calculate_z_level(index)
        return field_of_view, channel, z_level

    # ADDED BY ANDREW MARTENS MAY 2020
    # NOTE returns a stack in FORTRAN order! (good for ArrayFire)
    # NOTE modified from get_image below
    def get_image_stack(self, frame_number, field_of_view, z_level, height, width):
        """
        Attempts to get Image stack based on attributes alone.

        :type frame_number:  int
        :type field_of_view: int
        :type z_level:       int
        :type height:        int
        :type width:         int

        :rtype: Fortran-ordered Numpy array, or None
        # NOTE: dict mapping channel name to index will have to be inferred from somewhere else!!
        """
        # Use the chunk map to find the group of images, and read into a flat sequence of intensities
        image_group_number = self._calculate_image_group_number(
            frame_number, field_of_view, z_level
        )
        chunk = self._label_map.get_image_data_location(image_group_number)
        data = read_chunk(self._file_handle, chunk)

        # Convert the flat intensities to an "array"
        # NOTE: why not an ndarray?
        image_group_data = array.array("H", data)

        # NOTE What are the first 4 bytes? Timestamp? Is it bytes, or words, or ... ?
        image_data_start = 4

        # Based on the total length of the data, and the known width and height,
        # extrapolate the number of channels in the data.
        # NOTE why not use integer divide, // ?
        number_of_true_channels = (len(image_group_data) - image_data_start) // (
            height * width
        )
        # number_of_true_channels = int((len(image_group_data) - image_data_start) / (height * width))

        # Zeroth dimension is the channel
        # zeroth index -> first channel, 1st index -> 2nd channel etc.
        # This produces an array with channels as the first of 3 indices, which is inefficient for F-ordering.
        reshaped = np.reshape(
            image_group_data[image_data_start:],
            (number_of_true_channels, width, height),
            order="F",
        )

        # Fotran ordering works best when the last index is the slowest-changing
        # For image stacks, this typically means that 2D slices are stored in the 1st 2 dims,
        # and therefore a whole slice is in the 3rd dimension.
        # NOTE because the data are Fotran-ordered, they will appear Transposed
        # when viewing with matplotlib pyplot. However, it's more appropriate to Transpose later,
        # and to leave the data intact now.
        image_stack = np.empty(
            (width, height, number_of_true_channels), dtype=reshaped.dtype, order="F"
        )
        for c in range(number_of_true_channels):
            image_stack[..., c] = reshaped[c]

        # Skip images that are all zeros! This is important, since NIS Elements creates blank "gap" images if you
        # don't have the same number of images each cycle. We discovered this because we only took GFP images every
        # other cycle to reduce phototoxicity, but NIS Elements still allocated memory as if we were going to take
        # them every cycle.
        if np.any(image_stack):
            # NOTE: does this mean the first 8 bytes are the timestamp? 4 bytes? Why not parse out the timestamp first, then load the rest into the array? This was a goofy way of doing things by Rybarski...
            # NOTE do we care about the timestamp??
            # timestamp = struct.unpack("d", data[:8])[0]
            # return timestamp, image_stack
            return image_stack
        raise NoImageError

    def get_image(self, index):
        """
        Creates an Image object and adds its metadata, based on the index (which is simply the order in which the image was acquired). May return None if the ND2 contains
        multiple channels and not all were taken in each cycle (for example, if you take bright field images every minute, and GFP images every five minutes, there will be some
        indexes that do not contain an image. The reason for this is complicated, but suffice it to say that we hope to eliminate this possibility in future releases. For now,
        you'll need to check if your image is None if you're doing anything out of the ordinary.

        :type index:    int
        :rtype:    Image or None

        """
        field_of_view, channel, z_level = self.calculate_image_properties(index)
        channel_offset = index % len(self._metadata.channels)
        image_group_number = int(index / len(self._metadata.channels))
        frame_number = self._calculate_frame_number(
            image_group_number, field_of_view, z_level
        )
        try:
            timestamp, image = self._get_raw_image_data(
                image_group_number,
                channel_offset,
                self._metadata.height,
                self._metadata.width,
            )
        except NoImageError:
            return None
        else:
            image.add_params(
                index, timestamp, frame_number, field_of_view, channel, z_level
            )
            return image

    def get_image_by_attributes(
        self, frame_number, field_of_view, channel_name, z_level, height, width
    ):
        """
        Attempts to get Image based on attributes alone.

        :type frame_number:    int
        :type field_of_view:    int
        :type channel_name:    str
        :type z_level:    int
        :type height:    int
        :type width:    int

        :rtype: Image or None
        """
        image_group_number = self._calculate_image_group_number(
            frame_number, field_of_view, z_level
        )
        try:
            timestamp, raw_image_data = self._get_raw_image_data(
                image_group_number, self._channel_offset[channel_name], height, width
            )
            image = Image(raw_image_data)
            image.add_params(
                image_group_number,
                timestamp,
                frame_number,
                field_of_view,
                channel_name,
                z_level,
            )
        except (TypeError, NoImageError):
            return None
        else:
            return image

    def _calculate_field_of_view(self, index):
        """
        Determines what field of view was being imaged for a given image.

        :type index:    int
        :rtype:    int

        """
        images_per_cycle = len(self._metadata.z_levels) * len(self._metadata.channels)
        return int((index - (index % images_per_cycle)) / images_per_cycle) % len(
            self._metadata.fields_of_view
        )

    def _calculate_channel(self, index):
        """
        Determines what channel a particular image is.

        :type index:    int
        :rtype:    str

        """
        return self._metadata.channels[index % len(self._metadata.channels)]

    def _calculate_z_level(self, index):
        """
        Determines the plane in the z-axis a given image was taken in. In the future, this will be replaced with the actual offset in micrometers.

        :type index:    int
        :rtype:    int
        """
        return self._metadata.z_levels[
            int(
                (
                    (index - (index % len(self._metadata.channels)))
                    / len(self._metadata.channels)
                )
                % len(self._metadata.z_levels)
            )
        ]

    def _calculate_image_group_number(self, frame_number, fov, z_level):
        """
        Images are grouped together if they share the same time index, field of view, and z-level.

        :type frame_number: int
        :type fov: int
        :type z_level: int

        :rtype: int

        """
        return frame_number * len(self._metadata.fields_of_view) * len(
            self._metadata.z_levels
        ) + (fov * len(self._metadata.z_levels) + z_level)

    def _calculate_frame_number(self, image_group_number, field_of_view, z_level):
        """
        Images are in the same frame if they share the same group number and field of view and are taken sequentially.

        :type image_group_number:    int
        :type field_of_view:    int
        :type z_level:    int

        :rtype:    int

        """
        return (
            image_group_number
            - (field_of_view * len(self._metadata.z_levels) + z_level)
        ) / (len(self._metadata.fields_of_view) * len(self._metadata.z_levels))

    @property
    def _channel_offset(self):
        """
        Image data is interleaved for each image set. That is, if there are four images in a set, the first image
        will consist of pixels 1, 5, 9, etc, the second will be pixels 2, 6, 10, and so forth.

        :rtype: dict

        """
        return {channel: n for n, channel in enumerate(self._metadata.channels)}

    def _get_raw_image_data(self, image_group_number, channel_offset, height, width):
        """
        Reads the raw bytes and the timestamp of an image.

        :param image_group_number: groups are made of images with the same time index, field of view and z-level
        :type image_group_number: int
        :param channel_offset: the offset in the array where the bytes for this image are found
        :type channel_offset: int

        :rtype:    (int, Image)
        :raises:    NoImageError

        """
        chunk = self._label_map.get_image_data_location(image_group_number)
        data = read_chunk(self._file_handle, chunk)
        # print("data", data, "that was data")
        # All images in the same image group share the same timestamp! So if you have complicated image data,
        # your timestamps may not be entirely accurate. Practically speaking though, they'll only be off by a few
        # seconds unless you're doing something super weird.
        timestamp = struct.unpack("d", data[:8])[0]
        image_group_data = array.array("H", data)
        image_data_start = 4 + channel_offset

        # The images for the various channels are interleaved within the same array. For example, the second image
        # of a four image group will be composed of bytes 2, 6, 10, etc. If you understand why someone would design
        # a data structure that way, please send the author of this library a message.
        number_of_true_channels = int((len(image_group_data) - 4) / (height * width))
        image_data = np.reshape(
            image_group_data[image_data_start::number_of_true_channels], (height, width)
        )

        # Skip images that are all zeros! This is important, since NIS Elements creates blank "gap" images if you
        # don't have the same number of images each cycle. We discovered this because we only took GFP images every
        # other cycle to reduce phototoxicity, but NIS Elements still allocated memory as if we were going to take
        # them every cycle.
        if np.any(image_data):
            return timestamp, Image(image_data)
        raise NoImageError
