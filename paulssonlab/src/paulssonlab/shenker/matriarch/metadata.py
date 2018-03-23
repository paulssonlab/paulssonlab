import nd2reader
import PIL
import array
import zarr
from numcodecs import Blosc
import xmltodict
import xml
from util import recursive_map

DEFAULT_METADATA_COMPRESSOR = Blosc(
    cname="zstd", clevel=5, shuffle=Blosc.NOSHUFFLE, blocksize=0
)
ND2_METADATA_PARSED = [
    "image_metadata_sequence",
    "image_calibration",
    "image_attributes",
    "roi_metadata",
    "image_text_info",
    "image_metadata",
]
ND2_METADATA_XML = ["lut_data", "grabber_settings", "custom_data", "app_info"]
ND2_METADATA_ARRAYS = {
    "pfs_status": "i",
    "pfs_offset": "i",
    "x_data": "d",
    "y_data": "d",
    "z_data": "d",
    "camera_exposure_time": "d",
    "camera_temp": "d",
    "acquisition_times": "d",
    "acquisition_times_2": "d",
    "acquisition_frames": "i",
}  # TODO: have no idea about correct type for acquisition_frames
# NIKON_TIFF_METADATA_TAGS = [270, 65330, 65331, 65332, 65333]
NIKON_TIFF_METADATA_TAGS = [270, 65330, 65331, 65333]


def parse_nd2_file_metadata(nd2_file):
    return parse_nd2_metadata(nd2reader.ND2Reader(nd2_file))


def _stringify_dict_keys(d):
    return recursive_map(
        lambda s: s.decode("utf-8"), d, shortcircuit=bytes, ignore=True, keys=True
    )


def parse_nd2_metadata(nd2):
    label_map = nd2.parser._label_map
    raw_metadata = nd2.parser._raw_metadata
    metadata = {}
    metadata["time_source"] = _nd2_parse_chunk(nd2, b"CustomData|TimeSourceCache!")
    metadata["roi_rle_global"] = _nd2_parse_chunk(nd2, b"CustomData|RoiRleGlobal_v1!")
    metadata["acquisition_time"] = _nd2_parse_xml_chunk(
        nd2, b"CustomDataVar|AcqTimeV1_0!"
    )
    metadata["nd_control"] = _nd2_parse_xml_chunk(nd2, b"CustomDataVar|NDControlV1_0!")
    metadata["stream_data"] = _nd2_parse_xml_chunk(
        nd2, b"CustomDataVar|StreamDataV1_0!"
    )
    for label in ND2_METADATA_XML:
        location = getattr(label_map, label)
        if location:
            data = _nd2_parse_xml_chunk(nd2, location=location)
        else:
            data = None
        metadata[label] = data
    for label in ND2_METADATA_PARSED:
        metadata[label] = _stringify_dict_keys(getattr(raw_metadata, label))
    for label, dtype in ND2_METADATA_ARRAYS.items():
        metadata[label] = _nd2_parse_array(nd2, label, dtype)
    return metadata


def _nd2_parse_chunk(nd2, label=None, location=None):
    if location is None:
        if label is None:
            raise ValueError("need either label or location")
        location = nd2.parser._label_map._get_location(label)
    return nd2reader.common.read_chunk(nd2._fh, location)


def _nd2_parse_xml_chunk(nd2, label=None, location=None):
    data = _nd2_parse_chunk(nd2, label=label, location=location)
    if data is not None:
        return xmltodict.parse(data)
    else:
        return None


def _nd2_parse_array(nd2, label, dtype, compressor=DEFAULT_METADATA_COMPRESSOR):
    chunk_location = getattr(nd2.parser._label_map, label)
    raw_data = nd2reader.common.read_chunk(nd2._fh, chunk_location)
    if raw_data is None:
        return None
    raw_ary = array.array(dtype, raw_data)
    ary = zarr.array(raw_ary, compressor=compressor)
    return ary


def parse_nikon_tiff_file_metadata(tiff_file):
    with PIL.Image.open(tiff_file) as f:
        return parse_nikon_tiff_metadata(
            {tag: f.tag[tag][0] for tag in NIKON_TIFF_METADATA_TAGS if tag in f.tag}
        )


def parse_nikon_tiff_metadata(tags):
    metadata = {}
    for tag, data in tags.items():
        if tag == 270:
            if data:
                try:
                    parsed_data = xmltodict.parse(data)
                except xml.parsers.expat.ExpatError:
                    parsed_data = data
            else:
                parsed_data = ""
            metadata["image_description"] = parsed_data
        elif tag == 65330:
            # SEEMS TO STORE (with null bytes): CameraTemp1, Camera_ExposureTime1, PFS_OFFSET, PFS_STATUS
            label = _nikon_tiff_label(b"SLxImageTextInfo")
            idx = data.index(label) - 2
            md = nd2reader.common.read_metadata(data[idx:], 1)
            metadata["image_text_info"] = _stringify_dict_keys(md)
        elif tag == 65331:
            label = _nikon_tiff_label(b"SLxPictureMetadata")
            idx = data.index(label) - 2
            md = nd2reader.common.read_metadata(data[idx:], 1)
            metadata["image_metadata_sequence"] = _stringify_dict_keys(md)
        elif tag == 65332:
            # SEEMS TO STORE (UTF-16 encoded): AppInfo_V1_0, CustomDataV2_0, GrabberCameraSettingsV1_0, LUTDataV1_0
            label = _nikon_tiff_label(b"AppInfo_V1_0")
            idx = (
                data.index(label) + len(label) + 9
            )  # TODO: no idea what these bytes are
            # print(tag,label,idx,data[idx:idx+100])
            # from IPython import embed;embed()
            md = xmltodict.parse(data[idx:])
            metadata["app_info"] = md
        else:
            metadata[tag] = data
    return metadata


def _nikon_tiff_label(label):
    if type(label) != bytes:
        label = label.encode("utf-8")
    return b"\x00".join([label[i : i + 1] for i in range(len(label))])


def _nikon_tiff_field(label, data):
    blabel = _nikon_tiff_label(label)
    idx = data.index(blabel)
    return data[idx + len(blabel) : idx + 100]
    # subset = data[idx+len(blabel):]
    # idx2 = subset.index(b"\x00\x00")
    # return subset[:idx2]
