import itertools as it
from functools import compose, partial
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import zarr
from cytoolz import get_in

from paulssonlab.image_analysis.drift import find_feature_drift, get_drift_features
from paulssonlab.image_analysis.geometry import filter_rois, iter_roi_crops, shift_rois
from paulssonlab.image_analysis.segmentation.watershed import (
    segment as watershed_segment,
)
from paulssonlab.image_analysis.trench_detection import find_trenches
from paulssonlab.image_analysis.util import get_delayed


def crop_rois(img, rois):
    crops = {}
    # TODO: the islice is just for testing (we only deal with three trenches for FOV), otherwise every dask task takes a long time
    # for i, crop in it.islice(geometry.iter_roi_crops(img, rois), 100):
    for i, crop in iter_roi_crops(img, rois):
        crops[i] = crop
    return crops


# def segment_crops(crops):
#     masks = {}
#     for i, crop in crops.items():
#         masks[i] = segmentation.watershed.segment(crop)
#     return masks


# TODO: this is really boilerplatey, also we want finer task granularity than doing a whole FOV at once
# def measure_crops(label_images, intensity_images):
#     keys = label_images.keys() & intensity_images.keys()
#     return {k: measure_crop(label_images[k], intensity_images[k]) for k in keys}
def measure_crops(intensity_images):
    keys = intensity_images.keys()
    return {k: measure_crop(intensity_images[k]) for k in keys}


# def measure_crop(label_image, intensity_image):
# return pd.DataFrame(
#     skimage.measure.regionprops_table(
#         label_image,
#         intensity_image,
#         properties=(
#             "label",
#             "intensity_mean",
#         ),
#     )
# ).set_index("label")
def measure_crop(intensity_image):
    centerline = intensity_image[:, intensity_image.shape[1] // 2]
    return pd.Series(
        {
            # "p1": np.percentile(intensity_image, 1),
            # "p50": np.median(intensity_image),
            "p90": np.percentile(intensity_image, 90),
            # "p99": np.percentile(intensity_image, 99),
            # "mean": np.mean(intensity_image),
            # "centerline_mean": np.mean(centerline),
            # "centerline_median": np.median(centerline),
        },
        name="value",
    ).rename_axis(index="observable")


def measure_mask_crops(label_images):
    return {k: measure_mask_crop(v) for k, v in label_images.items()}


def measure_mask_crop(label_image):
    return pd.DataFrame(
        skimage.measure.regionprops_table(
            label_image,
            properties=(
                "label",
                "area",
                "axis_major_length",
                "axis_minor_length",
                "orientation",
                "centroid",
            ),
        )
    ).set_index("label")


def write_parquet(output_dir, measurements, position, t):
    df = pd.concat(
        {
            channel: pd.concat(channel_df, names=["roi_idx"])
            for channel, channel_df in measurements.items()
        },
        names=["channel"],
    ).reset_index()
    df["position"] = np.array(position).astype(np.uint16)
    df["t"] = np.array(t).astype(np.uint16)
    pq.write_to_dataset(
        pa.Table.from_pandas(df, preserve_index=False),
        Path(output_dir) / "measurements",
        partition_cols=["position", "t"],
        existing_data_behavior="delete_matching",
    )


def stack_dict(d, size=None):
    if size is None:
        size = max(d.keys()) + 1
    shape = next(iter(d.values())).shape
    null = np.full(shape, np.nan)
    return [d.get(idx, null) for idx in range(size)]


def _pad(ary, shape):
    return np.pad(
        ary,
        [(0, max(goal - current, 0)) for goal, current in zip(shape, ary.shape)],
        constant_values=np.nan,
    )


def write_zarr(filename, crops, t, max_t, channels):
    store = zarr.DirectoryStore(filename)  # DirectoryStoreV3(filename)
    if not filename.exists():
        num_rois = max(crops[channels[0]].keys()) + 1
        num_channels = len(channels)
        max_shape = np.max([crop.shape for crop in crops[channels[0]].values()], axis=0)
        shape = (num_rois, max_t, num_channels, *max_shape)
        chunks = (5, 1, num_channels, None, None)
        ary = zarr.open_array(
            store,
            mode="a",
            zarr_version=2,
            shape=shape,
            chunks=chunks,
            fill_value=np.nan,
        )
    else:
        ary = zarr.open_array(store, mode="a", zarr_version=2)
        max_shape = ary.shape[-2:]
    stack = np.array(
        [
            stack_dict(
                {
                    idx: _pad(crop.astype(np.float32), max_shape)
                    for idx, crop in crops[channel].items()
                },
                size=ary.shape[0],
            )
            for channel in channels
        ]
    ).swapaxes(0, 1)
    ary[:, t, ...] = stack


def process_fov(
    get_frame_func,
    position,
    ts,
    output_dir,
    segmentation_channel,
    measurement_channels,
    image_limits,
    find_trenches_kwargs={},
    dark=None,
    flats=None,
    delayed=True,
):
    delayed = get_delayed(delayed)
    channels = [
        segmentation_channel,
        *(set(measurement_channels) - set([segmentation_channel])),
    ]
    measurement_channels = measurement_channels
    rois = None
    shifts = {}
    write_tasks = []
    for prev_t, t in list(zip(it.chain([None], ts[:-1]), ts)):
        segmentation_img = delayed(get_frame_func)(position, segmentation_channel, t)
        if rois is None:
            rois = delayed(find_trenches)(
                segmentation_img, **{**dict(join_info=True), **find_trenches_kwargs}
            )
            shifts[t] = np.array([0, 0])
            initial_drift_features = delayed(get_drift_features)(
                segmentation_img, rois, shifts[t]
            )
        else:
            shifts[t] = delayed(find_feature_drift)(
                initial_drift_features,
                segmentation_img,
                rois,
                initial_shift2=shifts[prev_t],
            )
        shifted_rois = delayed(filter_rois)(
            delayed(shift_rois)(rois, shifts[t]), image_limits
        )
        crops = {}
        measurements = {}
        for channel in channels:
            if channel == segmentation_channel:
                crops[channel] = delayed(crop_rois)(segmentation_img, shifted_rois)
                # mask_crops = delayed(segment_crops)(crops[channel])
                # mask_measurements = delayed(measure_mask_crops)(mask_crops)
            else:
                img = delayed(get_frame_func)(position, channel, t)
                crops[channel] = delayed(crop_rois)(img, shifted_rois)
            if channel in measurement_channels:
                # measurements[channel] = delayed(measure_crops)(mask_crops, crops[channel])
                measurements[channel] = delayed(measure_crops)(crops[channel])
        metadata = dict(shifts=shifts)
        write_tasks.append(
            delayed(write_parquet)(output_dir, measurements, position, t)
        )
        # TODO
        max_t = 300
        write_tasks.append(
            delayed(write_zarr)(
                output_dir / f"crops_v={position}.zarr",
                crops,
                t,
                max_t,
                measurement_channels,
            )
        )
        # TODO: rois, metadata
    return write_tasks


def composite_for_segmentation(imgs):
    return np.sum([img / img.max() for img in imgs.values()], axis=0)


class DelayedStore:
    pass


class DelayedArrayStore(DelayedStore):
    pass


class DelayedTableStore(DelayedStore):
    pass


class CallbackQueue:
    def __init__(self, add_callback=None, before_callback=None, after_callback=None):
        self.add_callback = add_callback
        self.before_callback = before_callback
        self.after_callback = after_callback
        self.queue = []

    @classmethod
    def expand_conditions(conditions):
        return sum(
            [
                [(obj, k) for k in key] if isinstance(key, list) else [(obj, key)]
                for obj, key in conditions
            ],
            [],
        )

    def add(self, conditions, func, *args, **kwargs):
        for condition in conditions:
            if len(condition) != 2:
                raise ValueError(
                    f"expecting length 2 tuple (obj, key) instead of {condition}"
                )
        if args or kwargs:
            func = partial(func, *args, **kwargs)
        if self.add_callback:
            add_res = self.add_callback(conditions, func)
        self.queue.append((conditions, func, add_res))

    def run(self):
        while True:
            any_fired = False
            for idx in enumerate(self.queue):
                conditions, func, add_res = self.queue[idx]
                if all(key in obj for obj, key in self.expand_conditions(conditions)):
                    self._run(conditions, func, add_res)
                    any_fired = True

            if not any_fired:
                break
        return

    def _run(self, conditions, func, add_res):
        if self.before_callback:
            before_res = self.before_callback(conditions, func, add_res)
        else:
            before_res = None
        func()
        if self.after_callback:
            self.after_callback(conditions, func, add_res, before_res)


class Pipeline:
    DEFAULT_CONFIG = {
        "composite_func": composite_for_segmentation,
        "roi_detection_func": find_trenches,
        "track_drift": True,
        "segmentation_func": watershed_segment,
    }
    REQUIRED_CONFIG = ["segmentation_channels"]

    def __init__(self, output_dir, config=None, delayed=False):
        self.output_dir = Path(output_dir)
        if not config:
            config = {}
        config = {**self.DEFAULT_CONFIG, **config}
        self.config = config
        self.delayed = get_delayed(delayed)
        self.callbacks = CallbackQueue()

    def validate_config(self):
        for key in self.REQUIRED_CONFIG:
            if not get_in(self.config, *key):
                raise ValueError(f"missing required config key {key}")

    def add_callback(self, conditions, func, *args, **kwargs):
        return self.callbacks.add(conditions, self.delayed(func, *args, **kwargs))

    # we should pick a name that's better/more intuitive than handle_message
    def handle_message(self, msg):
        match msg:
            case {"type": "image", **info}:
                match info:
                    case {"image_type": "fish_barcode"}:
                        self.handle_fish_barcode(self, msg)
                    case _:
                        self.handle_image(self, msg)
            case {"type": "ome_metadata"}:
                print("got OME metadata")  # TODO
            case {"type": "nd2_metadata"}:
                print("got ND2 metadata")  # TODO: ??
            case {"type": "event", **info}:
                print("event", info)
            case {"type": "done"}:
                print("DONE")
            case _:
                # this exception should be caught, we don't want malformed messages to crash the self
                raise ValueError("cannot handle message", msg)

    def handle_image(self, msg):
        image = msg["image"]
        metadata = msg["metadata"]
        fov_num = metadata["fov_num"]
        t = metadata["t"]
        channel = metadata["channel"]
        # store raw image ("/raw")
        # preprocess image
        # store preprocessed whole frame image ("/frame")
        # if all segmentation channels available -> detect rois
        segmentation_frame_keys = [
            (fov_num, channel, t) for channel in self.config["segmentation_channels"]
        ]
        self.add_callback(
            [self.frames, segmentation_frame_keys], self.do_roi_detection, image
        )
        self.rois[(fov_num, t)] = self.delayed(
            compose(find_trenches, composite_for_segmentation), image
        )
        # USE FUNCS FROM CONFIG
        # if rois available -> crop
        # if crop segmentation available -> measure

    def do_roi_detection(self, image):
        pass

    def handle_fish_barcode(self, msg):
        pass  # TODO


# TODO: use a namedtuple (or typing.NamedTuple, or dataclass) for keys so that fields are named
def handle_image(self, msg):
    image = msg["image"]
    metadata = msg["metadata"]
    fov_num = metadata["fov_num"]
    t = metadata["t"]
    channel = metadata["channel"]
    raw_key = ("raw", fov_num, t, channel)
    # store raw image (in production, we won't do this, we will only store crops as we do below)
    self.array[raw_key] = image
    # TODO: we need a way to store per-frame metadata and write it to disk
    trenches_key = (
        "trenches",
        fov_num,
    )
    trenches = self.table.get(trenches_key)
    # check if we have done trench detection for this FOV
    if trenches is None and channel == trench_detection_channel:
        # if not, find trenches and save the resulting table
        trenches = self.delayed(new.image.find_trench_bboxes)(
            image, peak_func=trench_detection.peaks.find_peaks
        )
        self.table[trenches_key] = trenches
    # this list keeps track of all the raw frames that need to be cropped
    # frames for multiple channels will accumulate in this list until we get a frame for trench_detection_channel
    # if we have already processed such a frame, then keys_to_crop will contain only the current frame (raw_key)
    keys_to_crop = self.state.setdefault(("keys_to_crop", fov_num), [])
    keys_to_crop.append(raw_key)
    # we only can do further processing if we have already detected trenches for this FOV
    if trenches is not None:
        for raw_to_crop in keys_to_crop:
            crop_key = ("crops", *raw_to_crop[1:])
            # save trench crops for every frame in keys_to_crop
            self.array[crop_key] = self.delayed(crop_trenches)(
                self.array[raw_to_crop], trenches
            )
            segmentation_key = ("segmentation", fov_num, t, segmentation_channel)
            segmentation = self.array.get(segmentation_key)
            if segmentation is not None:
                if crop_key[-1] in measure_channels:
                    # if we have segmentation masks for this frame, we can immediately segment only this frame
                    keys_to_measure = [crop_key]
                else:
                    keys_to_measure = []
            else:
                # we don't have a segmentation mask yet, so we need to add to the keys_to_measure list
                keys_to_measure = self.state.setdefault(
                    ("keys_to_measure", fov_num, t), []
                )
                if crop_key[-1] in measure_channels:
                    # we want to measure this frame
                    keys_to_measure.append(crop_key)
                if crop_key[-1] == segmentation_channel:
                    # if this frame's channel is the segmentation channel, run segmentation
                    segmentation = self.delayed(segment_trenches)(self.array[crop_key])
                    self.array[segmentation_key] = segmentation
                    # once we have the segmentation mask, get measurements for the mask
                    self.table[
                        (
                            "mask_measurements",
                            *crop_key[1:],
                        )
                    ] = self.delayed(
                        measure_mask_crops
                    )(segmentation)
            segmentation = self.array.get(segmentation_key)
            # if we now have the segmentation mask, try measuring all frames in the keys_to_measure list
            if segmentation is not None:
                for crop_to_measure in keys_to_measure:
                    measurements_key = ("measurements", *crop_to_measure[1:])
                    self.table[measurements_key] = self.delayed(measure_crops)(
                        segmentation, self.array[crop_to_measure]
                    )
                self.state.pop(("keys_to_measure", fov_num, t), None)
        self.state.pop(("keys_to_crop", fov_num), None)
