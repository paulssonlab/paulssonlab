import numpy as np
import pandas as pd
from cytoolz import get_in
import cachetools
import nd2reader
import sys
from util import zip_dicts, multi_join, array_to_tuples, get_one, unzip_dicts
from metadata import parse_nd2_metadata
from diagnostics import expand_diagnostics_by_label


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


def _get_nd2_frame_list(sizes, channels):
    all_frames = [
        (filename, v, channel, t)
        for filename, file_sizes, file_channels in zip_dicts(sizes, channels)
        for v in range(file_sizes["v"])
        for channel in file_channels
        for t in range(file_sizes["t"])
    ]
    all_frames = pd.MultiIndex.from_tuples(all_frames)
    all_frames.names = ["filename", "position", "channel", "t"]
    # TODO: should we sort if we end up sorting below?
    all_frames, _ = all_frames.sortlevel()
    return all_frames


def get_nd2_frame_list(filenames):
    nd2s = {filename: get_nd2_reader(filename, memmap=False) for filename in filenames}
    sizes = {filename: nd2.sizes for filename, nd2 in nd2s.items()}
    metadata = {filename: parse_nd2_metadata(nd2) for filename, nd2 in nd2s.items()}
    parsed_metadata = {filename: nd2.metadata for filename, nd2 in nd2s.items()}
    channels = {filename: md["channels"] for filename, md in parsed_metadata.items()}
    all_frames = _get_nd2_frame_list(sizes, channels)
    positions = get_position_metadata(metadata)
    channels_to_idx = get_channels_to_indices(channels)
    all_frames = multi_join(all_frames, positions)
    all_frames = multi_join(all_frames, channels_to_idx)
    all_frames.sort_index(inplace=True)
    return all_frames, metadata, parsed_metadata


def concat_dataframes(dfs):
    df = pd.concat(dfs, sort=True)
    fields = get_one(dfs.keys())._fields
    df.index.names = fields + df.index.names[len(fields) :]
    return df


def concat_series(series):
    df = pd.concat(series, axis=1, sort=True).T
    df.index.names = get_one(series.keys())._fields
    return df


def trench_points_to_dataframe(trench_points):
    df = pd.concat(
        {
            trench_set: pd.DataFrame(
                {
                    "top": array_to_tuples(trench_set_points[0]),
                    "bottom": array_to_tuples(trench_set_points[1]),
                }
            )
            for trench_set, trench_set_points in trench_points.items()
        }
    )
    df.index.names = ["trench_set", "trench"]
    return df


def unzip_trench_info(trench_info):
    trench_points, trench_diag, trench_err = unzip_dicts(trench_info)
    trench_points = {
        idx: trench_points_to_dataframe(tp)
        for idx, tp in trench_points.items()
        if tp is not None
    }
    trench_points = concat_dataframes(trench_points)
    trench_diag = concat_series(trench_diag)
    trench_diag = expand_diagnostics_by_label(trench_diag)
    trench_diag.index.rename("trench_set", level=-1, inplace=True)
    return trench_points, trench_diag, trench_err
