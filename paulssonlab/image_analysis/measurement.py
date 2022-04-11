from .workflow import get_trench_bboxes
import pandas as pd


def _trench_diag_to_dataframe(trench_diag, sep="."):
    df = trench_diag.to_frame().T
    expanded_df = diagnostics.expand_diagnostics_by_label(df)
    expanded_df.index = expanded_df.index.droplevel(0)
    expanded_df.index.names = [*expanded_df.index.names[:-1], "trench_set"]
    return expanded_df


def _trench_info_to_dataframe(trench_info):
    trench_points, trench_diag, trench_err = trench_info
    if trench_err is not None:
        # TODO: write trench_err
        return None
    trench_diag = _trench_diag_to_dataframe(trench_info[1])
    # FROM: https://stackoverflow.com/questions/14744068/prepend-a-level-to-a-pandas-multiindex
    trench_diag = pd.concat([trench_diag], axis=1, keys=["diag"])
    trenches = pd.concat(
        [trench_points, util.multi_join(trench_info[0].index, trench_diag)], axis=1
    )
    return trenches


def _trenches_to_bboxes(trenches, image_limits):
    trench_bboxes = get_trench_bboxes(trenches, image_limits)
    if trench_bboxes is not None:
        trenches = pd.concat([trenches, trench_bboxes], axis=1)
    return trenches


find_trenches_diag = diagnostics.wrap_diagnostics(
    trench_detection.find_trenches, ignore_exceptions=True, pandas=True
)


def do_find_trenches(*key):
    frame = get_nd2_frame(*key)
    trench_info = find_trenches_diag(frame)
    return trench_info


def do_trenches_to_bboxes(trench_info, key=None, index_names=("filename", "position")):
    trenches = _trench_info_to_dataframe(trench_info)
    if trenches is None:
        return None
    if key is not None:
        trenches = pd.concat([trenches], names=index_names, keys=[key])
    trenches = _trenches_to_bboxes(trenches, image_limits=image_limits)
    return trenches


def do_get_trench_err(trench_info):
    trench_points, trench_diag, trench_err = trench_info
    if trench_err is None:
        return None
    if trench_points is not None:
        raise ValueError("expecting trench_points to be None")
    return trench_info


import pickle


def do_serialize_to_disk(
    data, filename, overwrite=True, skip_nones=True, format="pickle"
):
    if skip_nones:
        data = {k: v for k, v in data.items() if v is not None}
    if not overwrite and os.path.exists(filename):
        raise FileExistsError
    with open(filename, "wb") as f:
        if format == "arrow":
            buf = pa.serialize(data).to_buffer()
            f.write(buf)
        elif format == "pickle":
            pickle.dump(data, f)
    return data


def do_save_trenches(trenches, filename, overwrite=True):
    trenches = pd.concat(trenches)
    processing.write_dataframe_to_parquet(
        filename, trenches, merge=False, overwrite=overwrite
    )
    return trenches


def do_measure_and_write(trenches, frames, return_none=True, write=True, **kwargs):
    if trenches is None:
        return None
    trenches = filter_trenches(trenches)
    res = measure(trenches, frames, **kwargs)
    if write:
        processing.write_images_and_measurements(
            res,
            filename_func=filename_func,
            dataframe_format="parquet",
            write_images=True,
            write_measurements=True,
        )
    if return_none:
        return None
    else:
        return res


from dask.distributed import get_client


def split_task(size, func, data, *args, **kwargs):
    client = get_client()
    chunks = util.split_into(data, size)
    tasks = []
    for chunk in chunks:
        tasks.append(client.submit(func, chunk, *args, **kwargs))
    return client.gather(tasks)


def _measure(
    trenches,
    frames,
    measurement_func,
    segmentation_channel="MCHERRY",
    measure_channels=None,
    segmentation_func=trench_segmentation.segment,
    include_frame=True,
    frame_bits=8,
    frame_downsample=4,
    filename=None,
    position=None,
):
    frame_transformation = compose(
        processing.zarrify,
        partial(image.quantize, bits=frame_bits),
        partial(image.downsample, factor=frame_downsample),
    )
    trench_crops = processing._get_trench_crops(
        trenches,
        frames,
        include_frame=include_frame,
        frame_transformation=frame_transformation,
        filename=filename,
        position=position,
    )
    res = {}
    segmentation_masks = {}
    measurements = {}
    # segment
    for trench_set, crops_trench_channel_t in trench_crops.items():
        if trench_set == "_frame":
            continue
        for trench_idx, crops_channel_t in crops_trench_channel_t.items():
            for channel, crops_t in crops_channel_t.items():
                for t, crop in crops_t.items():
                    if measure_channels is not None and channel not in measure_channels:
                        continue
                    segmentation_key = (trench_set, trench_idx, segmentation_channel, t)
                    segmentation_mask = segmentation_masks.get(segmentation_key, None)
                    if segmentation_mask is None and segmentation_func is not None:
                        segmentation_mask = segmentation_func(
                            trench_crops[trench_set][trench_idx][segmentation_channel][
                                t
                            ]
                        )
                        segmentation_masks[segmentation_key] = segmentation_mask
                        # measure mask
                        if measurement_func is not None:
                            measurements[
                                ("mask", (trench_set, trench_idx, t))
                            ] = measurement_func(segmentation_mask, None)
                    # measure
                    if measurement_func is not None:
                        measurements[
                            (channel, (trench_set, trench_idx, t))
                        ] = measurement_func(segmentation_mask, crop)
    if measurement_func is not None:
        measurement_dfs = util.map_dict_levels(
            lambda k: (k[1], k[0], *k[2:]), measurements
        )
        for name, dfs in measurement_dfs.items():
            dfs = util.unflatten_dict(dfs)
            if isinstance(util.get_one(dfs, level=2), pd.Series):
                df = pd.concat(
                    {
                        channel: pd.concat(channel_dfs, axis=1).T
                        for channel, channel_dfs in dfs.items()
                    },
                    axis=1,
                )
            else:
                df = pd.concat(
                    {
                        channel: pd.concat(channel_dfs, axis=0)
                        for channel, channel_dfs in dfs.items()
                    },
                    axis=1,
                )
            df.index.names = ["trench_set", "trench", "t", *df.index.names[3:]]
            measurement_dfs[name] = df
        res["measurements"] = measurement_dfs
    images = dict(raw=trench_crops)
    if segmentation_func is not None:
        images["segmentation"] = util.unflatten_dict(segmentation_masks)
    res["images"] = images
    return res


measure = processing.iterate_over_groupby(["filename", "position"])(_measure)


def _trench_diag_to_dataframe(trench_diag, sep="."):
    df = trench_diag.to_frame().T
    expanded_df = diagnostics.expand_diagnostics_by_label(df)
    expanded_df.index = expanded_df.index.droplevel(0)
    expanded_df.index.names = [*expanded_df.index.names[:-1], "trench_set"]
    return expanded_df


def _trench_info_to_dataframe(trench_info):
    trench_points, trench_diag, trench_err = trench_info
    if trench_err is not None:
        # TODO: write trench_err
        return None
    trench_diag = _trench_diag_to_dataframe(trench_info[1])
    # FROM: https://stackoverflow.com/questions/14744068/prepend-a-level-to-a-pandas-multiindex
    trench_diag = pd.concat([trench_diag], axis=1, keys=["diag"])
    trenches = pd.concat(
        [trench_points, util.multi_join(trench_info[0].index, trench_diag)], axis=1
    )
    return trenches


def _trenches_to_bboxes(trenches, image_limits):
    trench_bboxes = get_trench_bboxes(trenches, image_limits)
    if trench_bboxes is not None:
        trenches = pd.concat([trenches, trench_bboxes], axis=1)
    return trenches


find_trenches_diag = diagnostics.wrap_diagnostics(
    trench_detection.find_trenches, ignore_exceptions=True, pandas=True
)


def do_find_trenches(*key):
    frame = get_nd2_frame(*key)
    trench_info = find_trenches_diag(frame)
    return trench_info


def do_trenches_to_bboxes(trench_info, key=None, index_names=("filename", "position")):
    trenches = _trench_info_to_dataframe(trench_info)
    if trenches is None:
        return None
    if key is not None:
        trenches = pd.concat([trenches], names=index_names, keys=[key])
    trenches = _trenches_to_bboxes(trenches, image_limits=image_limits)
    return trenches


def do_get_trench_err(trench_info):
    trench_points, trench_diag, trench_err = trench_info
    if trench_err is None:
        return None
    if trench_points is not None:
        raise ValueError("expecting trench_points to be None")
    return trench_info


import pickle


def do_serialize_to_disk(
    data, filename, overwrite=True, skip_nones=True, format="pickle"
):
    if skip_nones:
        data = {k: v for k, v in data.items() if v is not None}
    if not overwrite and os.path.exists(filename):
        raise FileExistsError
    with open(filename, "wb") as f:
        if format == "arrow":
            buf = pa.serialize(data).to_buffer()
            f.write(buf)
        elif format == "pickle":
            pickle.dump(data, f)
    return data


def do_save_trenches(trenches, filename, overwrite=True):
    trenches = pd.concat(trenches)
    processing.write_dataframe_to_parquet(
        filename, trenches, merge=False, overwrite=overwrite
    )
    return trenches


def do_measure_and_write(trenches, frames, return_none=True, write=True, **kwargs):
    if trenches is None:
        return None
    trenches = filter_trenches(trenches)
    res = measure(trenches, frames, **kwargs)
    if write:
        processing.write_images_and_measurements(
            res,
            filename_func=filename_func,
            dataframe_format="parquet",
            write_images=True,
            write_measurements=True,
        )
    if return_none:
        return None
    else:
        return res


def split_task(size, func, data, *args, **kwargs):
    client = get_client()
    chunks = util.split_into(data, size)
    tasks = []
    for chunk in chunks:
        tasks.append(client.submit(func, chunk, *args, **kwargs))
    return client.gather(tasks)


def run_analysis():
    save_trench_err_futures = {}


all_analysis_futures = {}
save_trenches_futures = {}
save_trench_err_futures = {}

all_trench_bboxes_futures = {}  # TODO: just for debugging

for filename, filename_frames in selected_frames.groupby("filename"):
    # analysis_futures = {}
    trench_bboxes_futures = {}
    trench_err_futures = {}
    for position, frames in filename_frames.groupby("position"):
        key = (filename, position)
        segmentation_channel = "CFP"  # TODO!!!
        frame_to_segment = frames.loc[
            IDX[:, :, [segmentation_channel], 0], :
        ]  # TODO: make pluggable
        trenches_future = client.submit(
            do_find_trenches, *frame_to_segment.index[0], priority=10
        )
        trench_err_futures[key] = client.submit(do_get_trench_err, trenches_future)
        trench_bboxes_future = client.submit(
            do_trenches_to_bboxes, trenches_future, (filename, position), priority=10
        )
        trench_bboxes_futures[key] = trench_bboxes_future
        all_trench_bboxes_futures[key] = trench_bboxes_future
        analysis_future = client.submit(
            split_task,
            5,
            do_measure_and_write,
            trench_bboxes_future,
            frames,
            measurement_func=_measurement_func,
            # measurement_func=None,
            # segmentation_func=None,
            segmentation_channel=segmentation_channel,  # TODO!!!
            return_none=True,
            write=True,
            priority=-10,
        )
        # analysis_futures[key] = analysis_future
        all_analysis_futures[key] = analysis_future
    # save trenches
    trenches_filename = filename_func(
        kind="trenches", extension="parquet", filename=filename
    )
    save_trenches_futures[filename] = client.submit(
        do_save_trenches,
        list(dict(sorted(trench_bboxes_futures.items())).values()),
        trenches_filename,
        priority=100,
    )
    trench_errs_filename = filename_func(
        kind="trench_errs", extension="pickle", filename=filename
    )
    save_trench_err_futures[filename] = client.submit(
        do_serialize_to_disk, trench_err_futures, trench_errs_filename, priority=100
    )
    # OPTIONALLY: stream analysis to master
