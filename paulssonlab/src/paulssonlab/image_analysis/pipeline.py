import itertools as it
from functools import reduce
from operator import and_, methodcaller
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import zarr
from cytoolz import compose, get_in
from skimage.filters import threshold_otsu
from skimage.measure import regionprops_table

from paulssonlab.image_analysis.delayed import (
    DelayedArrayStore,
    DelayedQueue,
    DelayedStore,
    DelayedTableStore,
)
from paulssonlab.image_analysis.drift import find_feature_drift, get_drift_features
from paulssonlab.image_analysis.geometry import (
    filter_and_shift_rois,
    get_image_limits,
    iter_roi_crops,
)
from paulssonlab.image_analysis.image import power_law_composite
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


def segment_crops(crops):
    masks = {}
    for i, crop in crops.items():
        masks[i] = watershed_segment(crop)
    return masks


# TODO: this is really boilerplatey, also we want finer task granularity than doing a whole FOV at once
def measure_crops(label_images, intensity_images):
    keys = label_images.keys() & intensity_images.keys()
    return {k: measure_crop(label_images[k], intensity_images[k]) for k in keys}


def measure_crop(label_image, intensity_image):
    return pd.DataFrame(
        regionprops_table(
            label_image,
            intensity_image,
            properties=(
                "label",
                "intensity_mean",
            ),
        )
    ).set_index("label")


def measure_mask_crops(label_images):
    return {k: measure_mask_crop(v) for k, v in label_images.items()}


def measure_mask_crop(label_image):
    return pd.DataFrame(
        regionprops_table(
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


def _common_keys(mappings):
    return list(reduce(and_, map(methodcaller("keys"), mappings)))


def measure_fish_crops(fish_images):
    channels = fish_images.keys()
    keys = _common_keys(fish_images.values())
    return {
        k: measure_fish_crop({ch: fish_images[ch][k] for ch in channels}) for k in keys
    }


def otsu_mean(x):
    return x[threshold_otsu(x) <= x].mean()


def measure_fish_crop(images):
    # centerline = intensity_image[:, intensity_image.shape[1] // 2]
    return pd.concat(
        {
            channel: pd.Series(
                {
                    "mean": np.mean(image),
                    "otsu_mean": otsu_mean(image),
                    # "p1": np.percentile(intensity_image, 1),
                    # "p50": np.median(intensity_image),
                    # "p90": np.percentile(intensity_image, 90),
                    # "p99": np.percentile(intensity_image, 99),
                    # "mean": np.mean(intensity_image),
                    # "centerline_mean": np.mean(centerline),
                    # "centerline_median": np.median(centerline),
                },
                name="value",
            ).rename_axis(index="observable")
            for channel, image in images.items()
        }
    )


class Pipeline:
    # TODO: move boilerplate here

    def handle_message(self, msg):
        self._handle_message(msg)
        self._queue.poll()


class DefaultPipeline(Pipeline):
    DEFAULT_CONFIG = {
        "composite_func": power_law_composite,
        "roi_detection_func": find_trenches,
        "correct_drift": True,
        "segmentation_func": watershed_segment,
        "segment": False,
        "measure": False,
        "measure_fish": True,
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
        self.first_t = {}
        self.last_t = {}
        self.initial_rois = DelayedTableStore(self._queue, output_dir / "rois")
        self.rois = DelayedTableStore(self._queue, output_dir / "rois")
        self.measurements = DelayedTableStore(self._queue, output_dir / "measurements")
        self.mask_measurements = DelayedTableStore(
            self._queue, output_dir / "mask_measurements"
        )
        self.raw_frames = DelayedArrayStore(self._queue, output_dir / "raw_frames")
        self.processed_frames = DelayedArrayStore(
            self._queue, output_dir / "processed_frames"
        )
        self.crops = DelayedArrayStore(self._queue, output_dir / "crops")
        self.segmentation_masks = DelayedArrayStore(
            self._queue, output_dir / "segmentation_masks"
        )
        self.initial_drift_features = DelayedStore(self._queue)
        self.image_limits = DelayedStore(self._queue)
        self.shifts = DelayedStore(self._queue)
        self.fish_raw_frames = DelayedArrayStore(
            self._queue, output_dir / "fish_raw_frames"
        )
        self.fish_processed_frames = DelayedArrayStore(
            self._queue, output_dir / "fish_processed_frames"
        )
        self.fish_crops = DelayedArrayStore(self._queue, output_dir / "fish_crops")
        self.fish_measurements = DelayedTableStore(
            self._queue, output_dir / "fish_measurements"
        )

    def _apply_default_config(self):
        config = {**self.DEFAULT_CONFIG, **self.config}
        config["trench_detection_channels"] = config.get(
            "trench_detection_channels"
        ) or config.get("segmentation_channels")
        if config["measure"] and not config["segment"]:
            raise ValueError("cannot enable measurement with segmentation disabled")
        self.config = config

    def delayed(self, func, *args, **kwargs):
        return self._queue.delayed(self._delayed(func), *args, **kwargs)

    def validate_config(self):
        for key in self.REQUIRED_CONFIG:
            if not get_in(self.config, *key):
                raise ValueError(f"missing required config key {key}")

    # we should pick a name that's better/more intuitive than handle_message
    def _handle_message(self, msg):
        match msg:
            case {"type": "image", **info}:
                match info:
                    case {"image_type": "fish_barcode"}:
                        self.handle_fish_barcode(msg)
                    case {"image_type": "science"}:
                        self.handle_image(msg)
                    case _:
                        raise ValueError(
                            "received image message with missing or unknown image_type"
                        )
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
        # print("&&&&&&")
        # print(self._queue._items)
        # print("&&&&&&")

    def handle_image(self, msg):
        # print("EEEE0", self._queue._items)
        image = msg["image"]
        metadata = msg["metadata"]
        # print()
        # print("*****************")
        # print()
        print("IMAGE", metadata)
        # print()
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
                self.delayed(
                    self.config["preprocess_func"],
                    raw_frame,
                    # channel=channel,
                    # **(self.config.get("preprocess_kwargs") or {}),
                ),
            )
            processed_frames = self.processed_frames
        else:
            # print("BYPASS")
            processed_frame = raw_frame
            processed_frames = self.raw_frames
        # if all segmentation channels available -> detect rois
        segmentation_frames = [
            processed_frames[(fov_num, seg_channel, t)]
            for seg_channel in self.config["segmentation_channels"]
        ]
        # get rois
        # TODO: segmentation_frame might be recomputed (undesirable) if multiple different
        # instances (created in different calls to _handle_image) are dependencies (args)
        # for different DelayedCallables that are themselves added to DelayedQueue
        # via DelayedStore.setdefault
        # segmentation_frame = self.segmentation_frames.setdefault(
        #     (fov_num, t),
        #     self.delayed(self.config["composite_func"], segmentation_frames),
        # )
        segmentation_frame = self.delayed(
            self.config["composite_func"], segmentation_frames
        )
        image_limits = self.image_limits.setdefault(
            fov_num,
            self.delayed(lambda img: get_image_limits(img.shape), processed_frame),
        )
        first_t = self.first_t.get(fov_num, None)
        last_t = self.last_t.get(fov_num, None)
        if last_t is None:
            last_t = t
        else:
            last_t = max(t, last_t)
        self.last_t[fov_num] = last_t
        if first_t is None or t == first_t:
            self.first_t[fov_num] = t
            initial_rois = self.initial_rois.setdefault(
                (fov_num, t),
                self.delayed(
                    self.config["roi_detection_func"],
                    segmentation_frame,
                    _keep_args=True,
                    **(self.config.get("roi_detection_kwargs") or {}),
                ),
            )
            shift = self.shifts.setdefault((fov_num, t), np.array([0, 0]))
            self.initial_drift_features.setdefault(
                fov_num,
                self.delayed(get_drift_features, segmentation_frame, initial_rois),
            )
        else:
            last_shift = self.shifts.get((fov_num, t - 1), None)
            initial_shift = self.shifts[(fov_num, first_t)]
            initial_rois = self.initial_rois[(fov_num, first_t)]
            shift = self.shifts.setdefault(
                (fov_num, t),
                self.delayed(
                    find_feature_drift,
                    self.initial_drift_features[fov_num],
                    segmentation_frame,
                    initial_rois,
                    initial_shift1=initial_shift,
                    initial_shift2=last_shift,
                ),
            )
        rois = self.rois.setdefault(
            (fov_num, t),
            self.delayed(
                filter_and_shift_rois,
                initial_rois,
                shift,
                image_limits,
            ),
        )
        # if rois available -> crop
        crops = self.crops.setdefault(
            (fov_num, channel, t), self.delayed(crop_rois, processed_frame, rois)
        )
        # if crop available -> segment crops
        if self.config["segment"]:
            seg_masks = self.segmentation_masks.setdefault(
                (fov_num, t),
                self.delayed(segment_crops, crops),
            )
        # if crop segmentation available -> measure
        if self.config["measure"]:
            self.mask_measurements.setdefault(
                (fov_num, t),
                self.delayed(measure_mask_crops, seg_masks),
            )
            self.measurements.setdefault(
                (fov_num, channel, t),
                self.delayed(measure_crops, seg_masks, crops),
            )

    def handle_fish_barcode(self, msg):
        image = msg["image"]
        metadata = msg["metadata"]
        print("FISH IMAGE", metadata)
        fov_num = metadata["fov_num"]
        t = metadata["t"]
        channel = metadata["channel"]
        # store raw image ("/raw")
        raw_frame = self.fish_raw_frames.setdefault((fov_num, channel, t), image)
        # preprocess image
        # store preprocessed whole frame image ("/frame")
        if self.config.get("preprocess_func"):
            processed_frame = self.fish_processed_frames.setdefault(
                (fov_num, channel, t),
                self.delayed(
                    self.config["preprocess_func"],
                    raw_frame,
                    # channel=channel,
                    # **(self.config.get("preprocess_kwargs") or {}),
                ),
            )
            processed_frames = self.fish_processed_frames
        else:
            processed_frame = raw_frame
            processed_frames = self.fish_raw_frames
        last_t = self.last_t.get(fov_num, None)
        if last_t is None:
            raise ValueError("received FISH image before ROI detection")
        rois = self.rois[(fov_num, last_t)]
        # if rois available -> crop
        self.fish_crops.setdefault(
            (fov_num, channel, t), self.delayed(crop_rois, processed_frame, rois)
        )
        if self.config["measure_fish"]:
            crops = {
                fish_channel: self.fish_crops[(fov_num, fish_channel, t)]
                for fish_channel in self.config["fish_measure_channels"]
            }
            self.fish_measurements.setdefault(
                (fov_num, t),
                self.delayed(measure_fish_crops, crops),
            )
