import nd2reader
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


def read_nd2_metadata(nd2):
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
