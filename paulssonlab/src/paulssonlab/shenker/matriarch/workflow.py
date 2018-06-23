import numpy as np
import pandas as pd
from cytoolz import get_in
import cachetools
import nd2reader
import sys


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
            filename: position_dataframe(
                [
                    p
                    for p in get_in(
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
                ]
            )
            for filename, md in metadata.items()
        }
    )
    positions.index.names = ["filename", "position"]
    return positions


def get_channels_to_indices(channels):
    df = pd.concat(
        {
            filename: pd.DataFrame(
                {"channel": channels, "channel_idx": range(len(channels))}
            )
            for filename, channels in channels.items()
        }
    )
    df = df.reset_index().set_index(["level_0", "channel"])[["channel_idx"]]
    df.index.names = ["filename", "channel"]
    return df


ND2READER_CACHE = cachetools.LFUCache(maxsize=48)
ND2_FRAME_CACHE = cachetools.LFUCache(maxsize=10**9, getsizeof=sys.getsizeof)


def _get_nd2_reader(filename, **kwargs):
    return nd2reader.ND2Reader(filename, **kwargs)


get_nd2_reader = cachetools.cached(cache=ND2READER_CACHE)(_get_nd2_reader)


def _get_nd2_frame(filename, position, channel, t, memmap=False):
    reader = get_nd2_reader(filename)
    # TODO: how slow is the channel lookup?
    channel_idx = reader.metadata["channels"].index(channel)
    # TODO: should I wrap in ndarray or keep as PIMS Frame?
    return np.array(reader.get_frame_2D(v=position, c=channel_idx, t=t, memmap=memmap))


get_nd2_frame_cached = cachetools.cached(cache=ND2_FRAME_CACHE)(_get_nd2_frame)


def get_nd2_frame(filename, position, channel, t, memmap=False, **kwargs):
    return get_nd2_frame_cached(filename, position, channel, t, memmap=memmap)
