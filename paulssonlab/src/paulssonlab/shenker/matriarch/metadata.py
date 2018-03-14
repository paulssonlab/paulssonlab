import nd2reader
import PIL
import zarr
from numcodecs import Blosc
import xmltodict

DEFAULT_METADATA_COMPRESSOR = Blosc(
    cname="zstd", clevel=5, shuffle=Blosc.NOSHUFFLE, blocksize=0
)
ND2_METADATA_PARSED = [
    "image_metadata_sequence",
    "image_calibration",
    "image_attributes",
    "lut_data",
    "grabber_settings",
    "custom_data",
    "app_info",
    "image_text_info",
]
ND2_METADATA_INT_ARRAYS = ["pfs_status", "pfs_offset"]
ND2_METADATA_DOUBLE_ARRAYS = [
    "x_data",
    "y_data",
    "z_data",
    "camera_exposure_time",
    "camera_temp",
    "acquisition_times",
    "acquisition_times_2",
    "acquisition_frames",
]
NIKON_TIFF_METADATA_TAGS = [65330, 65331, 65332, 65333]


def parse_nd2_file_metadata(nd2_file):
    return parse_nd2_metadata(nd2reader.ND2Reader(nd2_file))


def parse_nd2_metadata(nd2):
    label_map = nd2.parser._label_map
    raw_metadata = nd2.parser._raw_metadata
    metadata = {}
    metadata["time_source"] = _nd2_parse_chunk(nd2, b"CustomData|TimeSourceCache!")
    metadata["roi_rle_global"] = _nd2_parse_chunk(nd2, b"CustomData|RoiRleGlobal_v1!")
    metadata["acquisition_time"] = xmltodict.parse(
        _nd2_parse_chunk(nd2, b"CustomDataVar|AcqTimeV1_0!")
    )
    metadata["nd_control"] = xmltodict.parse(
        _nd2_parse_chunk(nd2, b"CustomDataVar|NDControlV1_0!")
    )
    metadata["stream_data"] = xmltodict.parse(
        _nd2_parse_chunk(nd2, b"CustomDataVar|StreamDataV1_0!")
    )
    for label in ND2_METADATA_PARSED:
        metadata[label] = getattr(raw_metadata, label)
    for label in ND2_METADATA_INT_ARRAYS:
        metadata[label] = _nd2_parse_array(nd2, label, "int")
    for label in ND2_METADATA_DOUBLE_ARRAYS:
        metadata[label] = _nd2_parse_array(nd2, label, "double")
    return metadata


def _nd2_parse_chunk(nd2, label):
    return nd2reader.common.read_chunk(
        nd2._fh, nd2.parser._label_map._get_location(label)
    )


def _nd2_parse_array(nd2, label, dtype, compressor=DEFAULT_METADATA_COMPRESSOR):
    raw_ary = nd2reader.common.read_array(
        nd2._fh, dtype, getattr(nd2.parser._label_map, label)
    )
    ary = zarr.array(raw_ary, compressor=compressor)
    return ary


def parse_nikon_tiff_file_metadata(tiff_file):
    with PIL.Image.open(tiff_file) as f:
        return parse_nikon_tiff_metadata(
            {tag: f.tag[tag][0] for tag in NIKON_TIFF_METADATA_TAGS}
        )


def parse_nikon_tiff_metadata(tags):
    metadata = {}
    for tag, data in tags.items():
        if tag == 65330:
            label = _nikon_tiff_label(b"SLxImageTextInfo")
            idx = data.index(label) - 2
            md = nd2reader.common.read_metadata(data[idx:], 1)
            metadata["image_text_info"] = md
        elif tag == 65331:
            label = _nikon_tiff_label("SLxPictureMetadata")
            idx = data.index(label) - 2
            md = nd2reader.common.read_metadata(data[idx:], 1)
            metadata["image_metadata_sequence"] = md
        elif tag == 65332:
            label = _nikon_tiff_label("AppInfo_V1_0")
            idx = (
                data.index(label) + len(label) + 9
            )  # TODO: no idea what these bytes are
            md = xmltodict.parse(data[idx:])
            metadata["app_info"] = md
        else:
            metadata[tag] = data
    return metadata


def _nikon_tiff_label(label):
    if type(label) != bytes:
        label = label.encode("utf-8")
    return b"\x00".join([label[i : i + 1] for i in range(len(label))])
