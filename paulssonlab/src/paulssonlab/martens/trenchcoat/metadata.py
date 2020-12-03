import nd2reader
import os
import tables
import numpy
import xmltodict
import xml.etree.ElementTree as ElementTree
import re
from pprint import pprint


def get_metadata(n):
    """
    Input a node in the HDF5 metadata section, return a dict
    with the following metadata:
    """
    metadata = {
        "channels": n.channels.read(),
        "fields_of_view": n.fields_of_view.read(),
        "frames": n.frames.read(),
        "width": n.width.read(),
        "height": n.height.read(),
        "pixel_microns": n.pixel_microns.read(),
        "unix_timestamp": n.unix_timestamp.read(),
        "z_levels": n.z_levels.read(),
    }

    return metadata


def make_fov_metadata_table_info_type():
    """
    Define the data types for a PyTables table.

    Table for looking up FOV, Frame, and returning the timestamp in
    seconds and X, Y, Z positions.
    """
    # Define the PyTables column types using a dictionary
    column_types = {
        "info_fov": tables.UInt16Col(),
        "info_frame": tables.UInt16Col(),
        "info_timestamp": tables.Float64Col(),
        "info_x": tables.Float64Col(),
        "info_y": tables.Float64Col(),
        "info_z": tables.Float64Col(),
        "info_pfs_offset": tables.Float64Col(),
        "info_pfs_status": tables.BoolCol(),
    }

    return type("FOV_Metadata", (tables.IsDescription,), column_types)


def dict_to_h5_metadata(dictionary, parent_node, h5file):
    """
    Recursively walk through a dictionary of metadata & copy over the
    metadata to an HDF5 file.
    """
    for k in dictionary.keys():
        name = k.decode("utf-8")
        # TODO: modify the name if it begins with: b (bool), ba (?), d (double), e (?), p (?), s (string), ui (unsigned integer), other? ->
        # basically, it seems the real name always begins with an uppercase letter, and is prefixed by
        # lowercase letters which have something to do with the data type?
        # This is similar, but not the same, as the types_xml, with the difference here that we can use numpy
        # to automatically do the type conversion for us (but not name conversion).

        # FIXME one of the variables is uiCon20(L --> is that an erroneous name? Seems weird!

        # Some of the variables are unnamed, so give it an explicit non-name
        if name == "":
            name = "no_name"

        # FIXME: doesn't work, perhaps because some names have different prefixes but then are the same.
        # else:
        ## Remove all letters until the first uppercase letter
        ## FIXME: performance boost if we compile the pattern just once & then pass it in?
        # pattern = re.compile("^[a-z]*(.*)")
        # try:
        # match = pattern.match(name)
        # name = match.group(1)
        ## If the pattern didn't match (e.g. was all lowercase)
        # if name == "":
        # name = match.group(0)
        # except:
        # print("Failed with name {}".format(name))

        elem = dictionary[k]

        if type(elem) == dict:
            new_group = h5file.create_group(parent_node, name)
            dict_to_h5_metadata(elem, new_group, h5file)

        elif type(elem) == list:
            new_group = h5file.create_group(parent_node, name)
            for number, obj in enumerate(elem):
                list_to_h5_metadata(obj, new_group, h5file, "Element_{}".format(number))

        else:
            new_node = h5file.create_array(parent_node, name, obj=numpy.array(elem))


def list_to_h5_metadata(elem, parent_node, h5file, name):
    """
    Recursively walk through a list of metadata & copy over the metadata to
    an HDF5 file.
    """
    if type(elem) == dict:
        new_group = h5file.create_group(parent_node, name)
        dict_to_h5_metadata(elem, new_group, h5file)

    elif type(elem) == list:
        new_group = h5file.create_group(parent_node, name)
        for number, obj in enumerate(elem):
            list_to_h5_metadata(obj, new_group, h5file, "Element_{}".format(number))

    else:
        new_node = h5file.create_array(parent_node, name, obj=numpy.array(elem))


def xml_to_h5_metadata_ascii(elem, parent_node, h5file, types_xml):
    """
    After re-converting the nested OrderedDicts (JSON-like) back into XML,
    iterate all the xml elements & copy over the metadata to an HDF5 file.

    NOTE: Unicode -> ASCII conversion to handle degree symbol and micron symbol
    (problems storing unicode symbols in pytables arrays?)
    """
    for child in elem:
        # Make a new group node
        if child.attrib["runtype"] == "CLxListVariant":
            new_group = h5file.create_group(parent_node, child.tag)
            xml_to_h5_metadata_ascii(child, new_group, h5file, types_xml)

        # Make a 1x1 carray data node
        else:
            # NOTE: converted to ascii
            # To decode, use decode('unicode-escape')
            new_node = h5file.create_array(
                parent_node,
                child.tag,
                obj=numpy.array(
                    child.attrib["value"].encode("ascii", "backslashreplace"),
                    dtype=types_xml[child.attrib["runtype"]],
                ),
            )


def copy_metadata(hdf5_dir, in_file, frames, fields_of_view):
    """
    Copy ND2 metadata into an HDF5 hierarchy.
    """
    reader = nd2reader.Nd2(in_file)
    (name, extension) = os.path.splitext(os.path.basename(in_file))

    out_file = os.path.join(hdf5_dir, "{}_metadata.h5".format(name))
    h5file = tables.open_file(out_file, mode="w")

    # Basic metadata
    copy_basic_metadata(h5file, reader, frames, fields_of_view)

    # Raw metadata
    copy_raw_metadata(h5file, reader)

    # FOV Metadata
    copy_fov_metadata(h5file, frames, fields_of_view, reader)

    h5file.close()


def copy_basic_metadata(h5file, reader, frames, fields_of_view):
    """
    Copy the "basic" metadata, which the ND2 library allows easy access for.
    """
    # Fields of view
    if not fields_of_view:
        fields_of_view = reader.fields_of_view

    h5file.create_array("/", "fields_of_view", obj=numpy.array(fields_of_view))

    # Z levels
    h5file.create_array("/", "z_levels", obj=numpy.array(reader.z_levels))

    # Frames
    if not frames:
        frames = reader.frames

    h5file.create_array("/", "frames", obj=numpy.array(frames))

    # Height
    h5file.create_array("/", "height", obj=numpy.array(reader.height))

    # Width
    h5file.create_array("/", "width", obj=numpy.array(reader.width))

    # Pixel Microns
    h5file.create_array("/", "pixel_microns", obj=numpy.array(reader.pixel_microns))

    # Date, in UNIX time (as a float)
    h5file.create_array("/", "unix_timestamp", obj=numpy.array(reader.date.timestamp()))

    # Channels
    h5file.create_array("/", "channels", obj=numpy.array(reader.channels))


def copy_raw_metadata(h5file, reader):
    """
    Copy the "raw" metadata, which are accessible but with more effort.
    """
    types_xml = {
        "CLxStringW": numpy.unicode,
        "lx_int32": numpy.int32,
        "lx_uint32": numpy.uint32,
        "double": numpy.float64,
        "bool": numpy.bool,
    }

    # Grabber Settings
    xml = xmltodict.unparse(reader._parser.raw_metadata.grabber_settings)
    tree = ElementTree.fromstring(xml)
    root_node = h5file.create_group(
        "/raw_metadata", "grabber_settings", createparents=True
    )
    xml_to_h5_metadata_ascii(tree, root_node, h5file, types_xml)

    # Image Attributes
    root_node = h5file.create_group(
        "/raw_metadata", "image_attributes", createparents=True
    )
    dict_to_h5_metadata(reader._parser.raw_metadata.image_attributes, root_node, h5file)

    # Image Calibration
    root_node = h5file.create_group(
        "/raw_metadata", "image_calibration", createparents=True
    )
    dict_to_h5_metadata(
        reader._parser.raw_metadata.image_calibration, root_node, h5file
    )

    # Image Metadata
    root_node = h5file.create_group(
        "/raw_metadata", "image_metadata", createparents=True
    )
    dict_to_h5_metadata(reader._parser.raw_metadata.image_metadata, root_node, h5file)

    # Image Metadata Sequence
    root_node = h5file.create_group(
        "/raw_metadata", "image_metadata_sequence", createparents=True
    )
    dict_to_h5_metadata(
        reader._parser.raw_metadata.image_metadata_sequence, root_node, h5file
    )

    # Image Text Info
    root_node = h5file.create_group(
        "/raw_metadata", "image_text_info", createparents=True
    )
    dict_to_h5_metadata(reader._parser.raw_metadata.image_text_info, root_node, h5file)

    # LUT Data
    xml = xmltodict.unparse(reader._parser.raw_metadata.lut_data)
    tree = ElementTree.fromstring(xml)
    root_node = h5file.create_group("/raw_metadata", "lut_data", createparents=True)
    xml_to_h5_metadata_ascii(tree, root_node, h5file, types_xml)

    # Custom Data
    xml = xmltodict.unparse(reader._parser.raw_metadata.custom_data)
    tree = ElementTree.fromstring(xml)
    root_node = h5file.create_group("/raw_metadata", "custom_data", createparents=True)
    xml_to_h5_metadata_ascii(tree, root_node, h5file, types_xml)


def copy_fov_metadata(h5file, frames, fields_of_view, reader):
    """
    Make a table of X, Y, Z, Timestamp, PFS information for each FOV.
    """
    # Check whether the user specified lists of frames & fields of view to process
    # (thereby excluding those not in the list)
    if not frames:
        frames = reader.frames

    if not fields_of_view:
        fields_of_view = reader.fields_of_view

    # FOV metadata: X, Y, Z, timestamp
    FOV_Frame_Time = make_fov_metadata_table_info_type()
    fov_metadata_table = h5file.create_table(
        "/", "fov_metadata", FOV_Frame_Time, "FOV Metadata"
    )
    fov_metadata_row = fov_metadata_table.row

    num_fov = len(reader.fields_of_view)
    for i, data in enumerate(
        zip(
            reader._parser.raw_metadata.x_data,
            reader._parser.raw_metadata.y_data,
            reader._parser.raw_metadata.z_data,
            reader._parser.raw_metadata.acquisition_times,
            reader._parser.raw_metadata.pfs_offset,
            reader._parser.raw_metadata.pfs_status,
        )
    ):

        # Convert linear index into FOV, frame
        fov = i % num_fov
        frame_number = int((i - fov) / num_fov)

        # Write the metadata into the table
        if (frame_number in frames) and (fov in fields_of_view):
            # NOTE index reader on FOV, b/c there might be gaps.
            # If no gaps, then reader.fields_of_view[fov] == fov.
            fov_metadata_row["info_fov"] = reader.fields_of_view[fov]
            fov_metadata_row["info_frame"] = frame_number
            fov_metadata_row["info_x"] = data[0]
            fov_metadata_row["info_y"] = data[1]
            fov_metadata_row["info_z"] = data[2]
            fov_metadata_row["info_timestamp"] = data[3]
            fov_metadata_row["info_pfs_offset"] = data[4]
            fov_metadata_row["info_pfs_status"] = data[5]

            fov_metadata_row.append()

    # Done populating the metadata table
    fov_metadata_table.flush()


def get_largest_extents(nd2_dict, metadata_key):
    """
    Input a dict of nd2 readers, and a numeric metadata key (fields_of_view,
    frames, z_levels...) Returns the smallest and largest values, across all of
    the files, as a tuple.

    Assume that values are non-negative integers. Assume values are
    always sorted from least to greatest.
    """
    from sys import maxsize

    smallest_value = maxsize
    largest_value = 0

    for _, reader in nd2_dict.items():
        values = getattr(reader, metadata_key)

        if values[0] < smallest_value:
            smallest_value = values[0]
        if values[-1] > largest_value:
            largest_value = values[-1]

    return (smallest_value, largest_value)


def compare_channel_names(nd2_list):
    """
    Input a list of nd2 file paths Returns false if not all files have
    identical channel name arrays in their metadata FIXME replace with the
    generic attributes function below?
    """
    iter_files = iter(nd2_list)
    zeroth_file_channels = nd2reader.Nd2(next(iter_files)).channels

    for f in iter_files:
        n = nd2reader.Nd2(f)
        if n.channels != zeroth_file_channels:
            return False

    return True


def metadata_attributes_equal(nd2_dict, attribute):
    """
    Input a dict of nd2 readers, and an attribute of interest Returns None
    if not all files have identical attributes of this type in their metadata
    Returns the attribute if they are identical.
    """
    iter_readers = iter(nd2_dict.items())
    zeroth_attribute = getattr(next(iter_readers)[1], attribute)

    for _, reader in iter_readers:
        if getattr(reader, attribute) != zeroth_attribute:
            return None

    return zeroth_attribute


def metadata_channels_equal(metadata_nodes):
    """
    Input list of nodes from the metadata node Returns false if not all
    nodes have identical channel name arrays in their metadata.
    """
    first_node = metadata_nodes.pop(0)
    first_channels = first_node.get_node("channels").read()

    for n in metadata_nodes:
        channels = n.get_node("channels").read()
        if not numpy.array_equal(first_channels, channels):
            return False

    return True


def main_print_metadata_function(in_file, sub_file_name):
    """
    Print the HDF5 metadata to the terminal. If no sub-file is specified, then print a list of all sub-files.
    This is useful, for example, if one wants to quickly see the names & order of the channels,
    or the image dimensions, etc.
    """
    h5file = tables.open_file(in_file, mode="r")

    # The user specified a sub-file
    if sub_file_name:
        node = h5file.get_node("/{}".format(sub_file_name))()
        metadata = get_metadata(node)
        pprint(metadata)
    # The user did not specify a sub-file
    else:
        for n in h5file.iter_nodes("/"):
            print(n._v_name)

    h5file.close()
