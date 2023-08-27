from pathlib import Path

import cachetools
import h5py
import nd2reader
import numpy as np
import pandas as pd
from cytoolz import get_in


def _get_position_list(md):
    return get_in(
        ["image_metadata", "SLxExperiment", "uLoopPars", "Points", ""], md
    ) or get_in(
        [
            "image_metadata",
            "SLxExperiment",
            "ppNextLevelEx",
            "",
            "uLoopPars",
            "Points",
            "",
        ],
        md,
    )


def get_position_metadata(metadata, grid_coords=True, reverse_grid="x"):
    def position_dataframe(d):
        df = pd.DataFrame.from_dict(d)
        df.rename(
            columns={
                "dPosName": "position_name",
                "dPosX": "x",
                "dPosY": "y",
                "dPosZ": "z",
                "dPFSOffset": "pfs_offset",
            },
            inplace=True,
        )
        df = df[["position_name", "x", "y", "z", "pfs_offset"]]
        if grid_coords:
            for coord in ("x", "y", "z"):
                coords = df[coord].unique()
                coords.sort()
                if coord in reverse_grid:
                    coords = coords[::-1]
                df[coord + "_idx"] = df[coord].map(
                    lambda c: np.where(coords == c)[0][0]
                )
        return df

    positions = pd.concat(
        {
            filename: position_dataframe(_get_position_list(md))
            for filename, md in metadata.items()
        }
    )
    positions.index.names = ["filename", "position"]
    return positions


ND2READER_CACHE = cachetools.LFUCache(maxsize=48)


class SplitFilename:
    def __init__(self, files):
        self.files = files

    def __str__(self):
        return self.files[0]

    def __repr__(self):
        return f"SplitFilename:{list(self.files)}"


def _get_nd2_reader(filename, **kwargs):
    if isinstance(filename, SplitFilename):
        from split_file_reader import SplitFileReader

        filename = SplitFileReader.open(filename.files, mode="rb")
    return nd2reader.ND2Reader(filename, **kwargs)


get_nd2_reader = cachetools.cached(cache=ND2READER_CACHE)(_get_nd2_reader)


def get_nd2_frame(filename, position, channel, t, dark=None, flat=None):
    reader = get_nd2_reader(filename)
    channel_idx = reader.metadata["channels"].index(channel)
    ary = reader.get_frame_2D(v=position, c=channel_idx, t=t)
    if dark is not None:
        ary = (
            ary - dark
        )  # can't subtract in place because img is uint16 and dark may be float64
    if flat is not None:
        ary = ary / flat
    return ary


def get_eaton_fish_frame(filename, v, channel, t, dark=None, flat=None):
    with h5py.File(Path(filename) / f"fov={v}_config={channel}_t={t}.hdf5") as f:
        ary = f["data"][()]
    if dark is not None:
        ary = (
            ary - dark
        )  # can't subtract in place because img is uint16 and dark may be float64
    if flat is not None:
        ary = ary / flat
    return ary


def get_filename_image_limits(metadata):
    return {
        filename: (
            (
                0,
                get_in(
                    ("image_attributes", "SLxImageAttributes", "uiWidth"), file_metadata
                )
                - 1,
            ),
            (
                0,
                get_in(
                    ("image_attributes", "SLxImageAttributes", "uiHeight"),
                    file_metadata,
                )
                - 1,
            ),
        )
        for filename, file_metadata in metadata.items()
    }


def points_dataframe(points, **kwargs):
    return pd.DataFrame({"x": points[:, 0], "y": points[:, 1]}, **kwargs)
