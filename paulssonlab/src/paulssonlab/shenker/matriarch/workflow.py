import numpy as np
import pandas as pd
import pyarrow as pa
import holoviews as hv
import streamz
from tornado import gen
from tornado.util import TimeoutError
from distributed.client import default_client
from cytoolz import get_in, valfilter, juxt, compose, partial
import cachetools
from numcodecs import Blosc
import nd2reader
import sys
from datetime import timedelta
from collections import defaultdict, namedtuple
from util import (
    zip_dicts,
    multi_join,
    array_to_tuples,
    get_one,
    unzip_dicts,
    iter_index,
    unzip_items,
    get_kwargs,
    kwcompose,
)
from metadata import parse_nd2_metadata
from geometry import get_image_limits, get_trench_bbox, bounding_box
from diagnostics import expand_diagnostics_by_label
from image import get_regionprops

IDX = pd.IndexSlice


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
ND2_FRAME_CACHE = cachetools.LFUCache(maxsize=10**8, getsizeof=sys.getsizeof)

# def _get_nd2_reader(filename, **kwargs):
def _get_nd2_reader(filename, memmap=False, **kwargs):
    return nd2reader.ND2Reader(filename, **kwargs)


# get_nd2_reader = _get_nd2_reader

get_nd2_reader = cachetools.cached(cache=ND2READER_CACHE)(_get_nd2_reader)
# get_nd2_reader = compose(lambda x: x.reopen(), _get_nd2_reader)


def _get_nd2_frame(filename, position, channel, t, memmap=False):
    reader = get_nd2_reader(filename)
    # TODO: how slow is the channel lookup?
    channel_idx = reader.metadata["channels"].index(channel)
    # ary = reader.get_frame_2D(v=position, c=channel_idx, t=t, memmap=memmap)
    ary = reader.get_frame_2D(v=position, c=channel_idx, t=t)
    # TODO: should I wrap in ndarray or keep as PIMS Frame?
    # return np.array(ary)
    return ary


get_nd2_frame = cachetools.cached(cache=ND2_FRAME_CACHE)(_get_nd2_frame)

get_nd2_frame_anyargs = kwcompose(
    get_nd2_frame, partial(get_kwargs, keys=["filename", "position", "channel", "t"])
)


def get_trench_image(
    trench_bboxes,
    filename,
    position,
    channel,
    t,
    trench_set,
    trench,
    *,
    get_frame_func=get_nd2_frame,
):
    frame = get_frame_func(filename, position, channel, t)
    trench_info = trench_bboxes.loc[
        IDX[filename, position, :, :, trench_set, trench], :
    ]
    # if len([s for s in trench_info.shape if s != 1]) != 1:
    #    raise Exception('trench_bboxes contains either more than one channel or more than one timepoint')
    uls = trench_info["upper_left"].values
    lrs = trench_info["lower_right"].values
    ul, lr = bounding_box(np.concatenate((uls, lrs)))
    return frame[ul[1] : lr[1] + 1, ul[0] : lr[0] + 1]


get_trench_set_image = partial(get_trench_image, trench=slice(None))


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
    all_frames["position_name"] = all_frames["position_name"].astype("category")
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


def points_dataframe(points, **kwargs):
    return pd.DataFrame({"x": points[:, 0], "y": points[:, 1]}, **kwargs)


def unzip_trench_info(trench_info):
    trench_points, trench_diag, trench_err = unzip_dicts(trench_info)
    trench_points = concat_dataframes(valfilter(lambda x: x is not None, trench_points))
    trench_diag = concat_series(trench_diag)
    trench_diag = expand_diagnostics_by_label(trench_diag)
    trench_diag.index.rename("trench_set", level=-1, inplace=True)
    return trench_points, trench_diag, trench_err


# def get_trench_thumbs(trenches, get_frame_func=get_nd2_frame):
#     trench_thumbs = {}
#     for idx, group in util.iter_index(trenches.groupby(['filename', 'position', 'channel', 't', 'trench_set'])):
#         frame = get_frame_func(**idx._asdict())
#         top_points = group['top'].values
#         bottom_points = group['bottom'].values
#         for trench_idx, _ in iter_index(group):
#             ul, lr = geometry.get_trench_bbox((top_points, bottom_points), trench_idx.trench, x_lim, y_lim)
#             trench_thumbs[trench_idx] = frame[ul[1]:lr[1],ul[0]:lr[0]]
#     return trench_thumbs


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


def _get_trench_bboxes(trenches, x_lim, y_lim, **kwargs):
    top_points = trenches["top"].values
    bottom_points = trenches["bottom"].values
    return np.hstack(
        [
            get_trench_bbox(
                (top_points, bottom_points), trench_idx, x_lim, y_lim, **kwargs
            )[:, np.newaxis]
            for trench_idx in range(len(top_points))
        ]
    )


def get_trench_bboxes(trenches, image_limits, **kwargs):
    def func(x):
        filename = x.index.get_level_values("filename")[0]
        x_lim, y_lim = image_limits[filename]
        upper_left, lower_right = _get_trench_bboxes(x, x_lim, y_lim, **kwargs)
        df = pd.concat(
            {
                "upper_left": points_dataframe(upper_left),
                "lower_right": points_dataframe(lower_right),
            },
            axis=1,
        )
        df.index.name = "trench"
        return df

    return trenches.groupby(
        ["filename", "position", "channel", "t", "trench_set"]
    ).apply(func)


def get_trench_stacks(
    trenches,
    frames,
    image_limits,
    get_frame_func=get_nd2_frame,
    transformation=np.array,
):
    trench_stacks = defaultdict(list)
    for framestack_idx, framestack_group in iter_index(
        trenches.groupby(["filename", "position", "channel"])
    ):
        # t_group_iterator = util.iter_index(framestack_group.groupby(['t', 'trench_set']))
        x_lim, y_lim = image_limits[framestack_idx.filename]
        current_trenches = next(iter(framestack_group.groupby("t")))[1]
        # optimize iter_index by not repeating steps in inner loops
        current_trenches_index = current_trenches.index.droplevel("t")
        Index = namedtuple("Index", current_trenches_index.names, rename=True)
        uls = current_trenches["upper_left"].values
        lrs = current_trenches["lower_right"].values
        for t in frames.loc[IDX[framestack_idx]].reset_index()["t"]:
            frame_idx = {"t": t, **framestack_idx._asdict()}
            frame = get_frame_func(**frame_idx)
            # for trench_idx, _ in iter_index(current_trenches_index):
            for idx in range(len(current_trenches_index)):
                trench_idx = Index(*current_trenches_index[idx])
                ul = uls[idx]
                lr = lrs[idx]
                trench_stacks[trench_idx].append(
                    frame[ul[1] : lr[1] + 1, ul[0] : lr[0] + 1]
                )
    if transformation is not None:
        trench_stacks = {k: transformation(v) for k, v in trench_stacks.items()}
    else:
        trench_stacks = dict(trench_stacks)  # convert back from a defaultdict
    return trench_stacks


def map_trenchwise(func, frame_stacks, trenches, channels=None):
    results = {}
    for trench_idx, _ in iter_index(trenches):
        if channels is None:
            results[trench_idx] = func(frame_stacks[trench_idx])
        else:
            results[trench_idx] = func(
                frame_stacks[trench_idx],
                *[
                    frame_stacks[trench_idx._replace(channel=channel)]
                    for channel in channels
                ],
            )
    return results


def map_stack(col_to_funcs, image_stack):
    ts = range(image_stack.shape[0])
    columns, funcs = unzip_items(col_to_funcs.items())
    func = juxt(*funcs)
    res = [func(image_stack[t]) for t in ts]
    d = {}
    for col, values in zip(columns, zip(*res)):
        if not isinstance(col, str):
            for i, sub_col in enumerate(col):
                d[sub_col] = [v[i] for v in values]
        else:
            d[col] = values
    df = pd.DataFrame(d, index=ts)
    df.index.name = "t"
    return df


def map_stack_over_labels(col_to_funcs, label_stack, intensity_stack, labels=None):
    dfs = {
        t: map_frame_over_labels(
            col_to_funcs, label_stack[t], intensity_stack[t], labels=labels
        )
        for t in range(label_stack.shape[0])
    }
    df = pd.concat(dfs)
    df.index.set_names("t", level=0, inplace=True)
    return df


# TODO
# slightly faster alternative to, e.g.,
# pd.DataFrame({'label': label_image.ravel(), 'value': intensity_image.ravel()}).groupby('label').agg(['mean', 'min', 'max'])
# OR numpy_indexed.group_by(l0.ravel(), i0.ravel(), reduction=np.mean)
def map_frame_over_labels(col_to_funcs, label_image, intensity_image, labels=None):
    if labels is None:
        labels = range(0, np.max(np.asarray(label_image)) + 1)
    columns, funcs = unzip_items(col_to_funcs.items())
    func = juxt(*funcs)
    res = [func(intensity_image[label_image == label]) for label in labels]
    d = {}
    for col, values in zip(columns, zip(*res)):
        if not isinstance(col, str):
            for i, sub_col in enumerate(col):
                d[sub_col] = [v[i] for v in values]
        else:
            d[col] = values
    df = pd.DataFrame(d, index=labels)
    df.index.name = "label"
    return df


def map_frame(col_to_funcs, image):
    columns, funcs = unzip_items(col_to_funcs.items())
    func = juxt(*funcs)
    res = func(image)
    d = {}
    for col, value in zip(columns, res):
        if not isinstance(col, str):
            for i, sub_col in enumerate(col):
                d[sub_col] = value[i]
        else:
            d[col] = value
    return pd.Series(d)


def _analyze_trench(
    trench_idx,
    frames_to_analyze,
    trench_images,
    trenchwise_funcs=None,
    labelwise_funcs=None,
    regionprops=False,
    segment_func=None,
):
    segmentation_channel = trench_idx.channel
    readout_channels = set(frames_to_analyze.index.get_level_values("channel"))
    if trenchwise_funcs:
        trenchwise_df = pd.concat(
            {
                channel: map_frame(trenchwise_funcs, trench_images[channel])
                for channel in readout_channels
            },
            axis=0,
        )
    else:
        trenchwise_df = None
    if labelwise_funcs or regionprops:
        label_trench_image = segment_func(trench_images[segmentation_channel])
        label_dfs = {}
        for channel in readout_channels:
            label_channel_dfs = {}
            if labelwise_funcs:
                label_channel_dfs["labelwise"] = map_frame_over_labels(
                    labelwise_funcs, label_trench_image, trench_images[channel]
                )
            if regionprops:
                regionprops_df = get_regionprops(
                    label_trench_image, trench_images[channel]
                )
                if regionprops_df is not None:
                    label_channel_dfs["regionprops"] = regionprops_df
            label_dfs[channel] = pd.concat(label_channel_dfs, axis=1)
        labelwise_df = pd.concat(label_dfs, axis=1)
        labelwise_df.index.name = "label"
    else:
        labelwise_df = None
    return trenchwise_df, labelwise_df


def analyze_trenches(
    trenches,
    frames_to_analyze,
    framewise_funcs=None,
    trenchwise_funcs=None,
    labelwise_funcs=None,
    regionprops=False,
    segment_func=None,
    get_frame_func=get_nd2_frame,
):
    frame_t_idx = tuple([frames_to_analyze.index[0][i] for i in (0, 1, 3)])
    readout_channels = set(frames_to_analyze.index.get_level_values("channel"))
    channels = {trenches.index.get_level_values("channel")[0], *readout_channels}
    uls = trenches["upper_left"].values
    lrs = trenches["lower_right"].values
    images = {
        channel: get_frame_func(
            filename=frame_t_idx[0],
            position=frame_t_idx[1],
            channel=channel,
            t=frame_t_idx[2],
        )
        for channel in channels
    }
    if framewise_funcs:
        framewise_df = pd.concat(
            {
                channel: map_frame(framewise_funcs, images[channel])
                for channel in readout_channels
            },
            axis=0,
        )
    else:
        framewise_df = None
    trenchwise_dfs = {}
    labelwise_dfs = {}
    for trench_idx, idx in iter_index(trenches.index):
        ul = uls[idx.Index]
        lr = lrs[idx.Index]
        trench_images = {
            channel: images[channel][ul[1] : lr[1] + 1, ul[0] : lr[0] + 1]
            for channel in images.keys()
        }
        trench_trenchwise_df, trench_labelwise_df = _analyze_trench(
            trench_idx,
            frames_to_analyze,
            trench_images,
            trenchwise_funcs=trenchwise_funcs,
            labelwise_funcs=labelwise_funcs,
            regionprops=regionprops,
            segment_func=segment_func,
        )
        result_idx = (*frame_t_idx, trench_idx.trench_set, trench_idx.trench)
        trenchwise_dfs[result_idx] = trench_trenchwise_df
        labelwise_dfs[result_idx] = trench_labelwise_df
    framewise_df = pd.concat({frame_t_idx: framewise_df}, axis=1).T
    trenchwise_df = pd.concat(trenchwise_dfs, axis=1).T
    labelwise_df = pd.concat(labelwise_dfs, axis=0)
    framewise_df.index.names = [
        "filename",
        "position",
        "t",
        *framewise_df.index.names[3:],
    ]
    for df in (trenchwise_df, labelwise_df):
        df.index.names = [
            "filename",
            "position",
            "t",
            "trench_set",
            "trench",
            *df.index.names[5:],
        ]
    return framewise_df, trenchwise_df, labelwise_df


def analyze_frames_and_trenches_iter(selected_trenches, frames_to_analyze, func):
    frames_t = frames_to_analyze.groupby(["filename", "position"])
    res = []
    for frame_idx, trenches in iter_index(
        selected_trenches.groupby(["filename", "position"])
    ):
        if not (frame_idx.filename, frame_idx.position) in frames_t.groups:
            continue
        fp_frames = frames_t.get_group((frame_idx.filename, frame_idx.position))
        for t, frames in fp_frames.groupby("t"):
            yield func(trenches, frames)
    return res


# TODO
def analyze_frames_and_trenches(selected_trenches, frames_to_analyze, func):
    frames_t = frames_to_analyze.groupby(["filename", "position"])
    res = []
    for frame_idx, trenches in iter_index(
        selected_trenches.groupby(["filename", "position"])
    ):
        if not (frame_idx.filename, frame_idx.position) in frames_t.groups:
            continue
        fp_frames = frames_t.get_group((frame_idx.filename, frame_idx.position))
        for t, frames in fp_frames.groupby("t"):
            res.append(func(trenches, frames))
    return res


def concat_unzip_dataframes(res):
    return [
        pd.concat(filter(lambda x: x is not None, dfs), axis=0) for dfs in zip(*res)
    ]


def sink_to_arrow(batches, sinks, writers, output_func=None):
    if output_func is None:
        output_func = lambda i: pa.BufferOutputStream()
    for i, batch in enumerate(batches):
        if i not in writers:
            sinks[i] = output_func(i)
            writers[i] = pa.RecordBatchStreamWriter(sinks[i], batch.schema)
        writers[i].write_batch(batch)


@streamz.Stream.register_api()
class gather_and_cancel(streamz.Stream):
    def __init__(
        self,
        upstream,
        stream_name=None,
        client=None,
        gather=True,
        cancel=True,
        timeout=None,
        timeout_func=None,
        success_func=None,
        loop=None,
    ):
        if client is None:
            client = default_client()
        if loop is None:
            loop = client.loop
        self.client = client
        self.gather = gather
        self.cancel = cancel
        if not isinstance(timeout, timedelta):
            timeout = timedelta(seconds=timeout)
        self.timeout = timeout
        self.timeout_func = timeout_func
        self.success_func = success_func
        streamz.Stream.__init__(self, upstream, stream_name=stream_name, loop=loop)

    @gen.coroutine
    def update(self, x, who=None):
        if self.gather:
            if self.timeout is not None:
                try:
                    result = yield gen.with_timeout(
                        self.timeout, self.client.gather(x, asynchronous=True)
                    )
                except TimeoutError:
                    if self.timeout_func is not None:
                        self.timeout_func(x)
                    # TODO: what happens if we get rid of this?
                    raise gen.Return(None)
            else:
                result = yield self.client.gather(x, asynchronous=True)
        else:
            result = x
        if self.cancel:
            yield self.client.cancel(x, asynchronous=True)
        if self.success_func is not None:
            self.success_func(x)
        result2 = yield self._emit(result)
        raise gen.Return(result2)


@streamz.Stream.register_api()
class with_timeout(streamz.Stream):
    def __init__(
        self,
        upstream,
        stream_name=None,
        timeout=None,
        retries=1,
        failure_func=None,
        loop=None,
    ):
        if not isinstance(timeout, timedelta):
            timeout = timedelta(seconds=timeout)
        self.timeout = timeout
        self.retries = retries
        self.failure_func = failure_func
        streamz.Stream.__init__(self, upstream, stream_name=stream_name, loop=loop)

    @gen.coroutine
    def update(self, x, who=None):
        future = self._emit(x)
        if self.timeout is not None:
            future = gen.with_timeout(self.timeout, future)
        got_result = False
        for attempt in range(self.retries):
            try:
                result = yield future
                got_result = True
                break
            except TimeoutError:
                continue
        if got_result:
            raise gen.Return(result)
        else:
            self.failure_func(x)
            raise gen.Return(None)


async def gather_stream(source, as_completed):
    async for future in as_completed:
        await source.emit(future)


### TEST
import time, gc


def _analyze_trench_dummy(
    trench_idx,
    frames_to_analyze,
    trench_images,
    trenchwise_funcs=None,
    labelwise_funcs=None,
    regionprops=False,
    segment_func=None,
):
    segmentation_channel = trench_idx.channel
    readout_channels = set(frames_to_analyze.index.get_level_values("channel"))
    if trenchwise_funcs:
        trenchwise_df = {
            channel: map_frame(trenchwise_funcs, trench_images[channel])
            for channel in readout_channels
        }
    else:
        trenchwise_df = None
    # time.sleep(0.01)
    # return
    if labelwise_funcs or regionprops:
        label_trench_image = segment_func(trench_images[segmentation_channel])
        time.sleep(0.01)
        return
        for channel in readout_channels:
            if labelwise_funcs:
                _ = map_frame_over_labels_dummy(
                    labelwise_funcs, label_trench_image, trench_images[channel]
                )
            if regionprops:
                regionprops_df = get_regionprops(
                    label_trench_image, trench_images[channel]
                )
    return None


def analyze_trenches_dummy(
    trenches,
    frames_to_analyze,
    framewise_funcs=None,
    trenchwise_funcs=None,
    labelwise_funcs=None,
    regionprops=False,
    segment_func=None,
):
    frame_t_idx = tuple([frames_to_analyze.index[0][i] for i in (0, 1, 3)])
    readout_channels = set(frames_to_analyze.index.get_level_values("channel"))
    channels = {trenches.index.get_level_values("channel")[0], *readout_channels}
    uls = trenches["upper_left"].values
    lrs = trenches["lower_right"].values
    # TEST1
    # time.sleep(30)
    # return None
    images = {}
    for channel in channels:
        # nd2 = nd2reader.ND2Reader(frame_t_idx[0])
        # c = nd2.metadata['channels'].index(channel)
        # img = nd2.get_frame_2D(v=frame_t_idx[1],
        #                       c=c,
        #                       t=frame_t_idx[2])
        img = np.load("/home/jqs1/projects/matriarch/temp/test.npy")
        images[channel] = img
    #     images = {channel: get_frame_func(filename=frame_t_idx[0],
    #                                       position=frame_t_idx[1],
    #                                       channel=channel,
    #                                       t=frame_t_idx[2]) for channel in channels}
    # TEST2
    # time.sleep(5)
    # time.sleep(25)
    # return None
    for trench_idx, idx in iter_index(trenches.index):
        ul = uls[idx.Index]
        lr = lrs[idx.Index]
        trench_images = {
            channel: images[channel][ul[1] : lr[1] + 1, ul[0] : lr[0] + 1]
            for channel in images.keys()
        }
        # TEST3
        _analyze_trench_dummy(
            trench_idx,
            frames_to_analyze,
            trench_images,
            trenchwise_funcs=trenchwise_funcs,
            labelwise_funcs=labelwise_funcs,
            regionprops=regionprops,
            segment_func=segment_func,
        )
    gc.collect()
    return None


def map_frame_over_labels_dummy(
    col_to_funcs, label_image, intensity_image, labels=None
):
    if labels is None:
        labels = range(0, np.max(np.asarray(label_image)) + 1)
    columns, funcs = unzip_items(col_to_funcs.items())
    func = juxt(*funcs)
    res = [func(intensity_image[label_image == label]) for label in labels]
    d = {}
    for col, values in zip(columns, zip(*res)):
        if not isinstance(col, str):
            for i, sub_col in enumerate(col):
                d[sub_col] = [v[i] for v in values]
        else:
            d[col] = values
    # df = pd.DataFrame(d, index=labels)
    # df.index.name = 'label'
    return d


def select_dataframe(df, params, **kwargs):
    params = {**params, **kwargs}
    idx = tuple(params.get(column, slice(None)) for column in df.index.names)
    # TODO: omitting any trailing slice(None) yields massive speedup
    while len(idx) and idx[-1] == slice(None):
        idx = idx[:-1]
    if any(isinstance(obj, slice) for obj in idx):
        # TODO: slices guarantee that we won't drop levels in .loc
        result = df.loc[idx, :]
    else:
        # TODO: .loc would drop levels
        result = df.xs(idx, drop_level=False)
    if isinstance(result, pd.Series):
        result = result.to_frame(name=0).T
    return result
