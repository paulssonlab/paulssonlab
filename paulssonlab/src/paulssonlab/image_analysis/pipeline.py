import itertools as it
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import zarr
from cytoolz import compose, get_in

from paulssonlab.image_analysis.delayed import (
    DelayedArrayStore,
    DelayedQueue,
    DelayedTableStore,
)
from paulssonlab.image_analysis.drift import find_feature_drift, get_drift_features
from paulssonlab.image_analysis.geometry import filter_rois, iter_roi_crops, shift_rois
from paulssonlab.image_analysis.image import mean_composite
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


class Pipeline:
    pass  # TODO: move boilerplate here


class DefaultPipeline(Pipeline):
    DEFAULT_CONFIG = {
        "composite_func": mean_composite,
        "roi_detection_func": find_trenches,
        "track_drift": True,
        "segmentation_func": watershed_segment,
    }
    REQUIRED_CONFIG = ["segmentation_channels"]

    def __init__(self, output_dir, config=None, delayed=False):
        self.output_dir = Path(output_dir)
        if not config:
            config = {}
        self.config = config
        self._apply_default_config()
        self._delayed = get_delayed(delayed)
        self._queue = DelayedQueue()
        self.rois = DelayedTableStore(self._queue, output_dir / "rois")
        self.measurements = DelayedTableStore(self._queue, output_dir / "measurements")
        self.raw_frames = DelayedArrayStore(self._queue, output_dir / "raw")
        self.processed_frames = DelayedArrayStore(self._queue, output_dir / "processed")
        self.crops = DelayedArrayStore(self._queue, output_dir / "crops")
        self.segmentation_masks = DelayedArrayStore(
            self._queue, output_dir / "segmentation_masks"
        )
        self.fish_crops = DelayedArrayStore(self._queue, output_dir / "fish_crops")

    def _apply_default_config(self):
        config = {**self.DEFAULT_CONFIG, **self.config}
        config["trench_detection_channels"] = config.get(
            "trench_detection_channels"
        ) or config.get("segmentation_channels")
        self.config = config

    def delayed(self, func, *args, **kwargs):
        return self._queue.delayed(self._delayed(func), *args, **kwargs)

    def validate_config(self):
        for key in self.REQUIRED_CONFIG:
            if not get_in(self.config, *key):
                raise ValueError(f"missing required config key {key}")

    # we should pick a name that's better/more intuitive than handle_message
    def handle_message(self, msg):
        match msg:
            case {"type": "image", **info}:
                match info:
                    case {"image_type": "fish_barcode"}:
                        self.handle_fish_barcode(msg)
                    case _:
                        self.handle_image(msg)
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
        print("IMAGE", metadata)
        fov_num = metadata["fov_num"]
        t = metadata["t"]
        channel = metadata["channel"]
        # store raw image ("/raw")
        raw_frame = self.raw_frames.setdefault((fov_num, channel, t), image)
        # preprocess image
        # store preprocessed whole frame image ("/frame")
        if self.config.get("preprocess_func"):
            processed_frame = self.processed_frames.setdefault(
                (fov_num, channel, t),
                self.delayed(self.config["preprocess_func"], image),
            )
        else:
            processed_frame = raw_frame
        # if all segmentation channels available -> detect rois
        segmentation_frame_keys = [
            (fov_num, seg_channel, t)
            for seg_channel in self.config["segmentation_channels"]
        ]
        # get rois
        roi_detection_func = compose(
            self.config["roi_detection_func"], self.config["composite_func"]
        )
        segmentation_frames = [self.raw_frames[k] for k in segmentation_frame_keys]
        rois = self.rois.setdefault(
            (fov_num, t), self.delayed(roi_detection_func, segmentation_frames)
        )
        # if rois available -> crop
        crops = self.crops.setdefault(
            (fov_num, channel, t), self.delayed(crop_rois, rois, processed_frame)
        )
        # # if crop available -> segment crops
        # seg_masks = self.segmentation_masks.setdefault(
        #     (fov_num, t),
        #     self.delayed(segment_crops, rois, crops,
        # )
        # # if crop segmentation available -> measure
        ### TODO: segmentationless measurement
        self.measurements.setdefault(
            (fov_num, t, channel),
            self.delayed(measure_crops, crops),
            # seg_masks,
        )

    def handle_fish_barcode(self, msg):
        pass  # TODO
