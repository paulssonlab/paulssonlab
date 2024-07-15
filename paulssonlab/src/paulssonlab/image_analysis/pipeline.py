from functools import partial, reduce
from operator import and_, methodcaller
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
from cytoolz import compose, get_in
from skimage.filters import threshold_otsu
from skimage.measure import regionprops_table
from skimage.registration import phase_cross_correlation

from paulssonlab.image_analysis.delayed import (
    DelayedBatchedZarrStore,
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


def combine_crops(func, crop_sets, *args, **kwargs):
    keys = _common_keys(crop_sets)
    return {
        k: func([crop_set[k] for crop_set in crop_sets], *args, **kwargs) for k in keys
    }


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


def measure_fish_crop(images):
    # TODO: can surely do this faster/more elegantly
    composite = power_law_composite(list(images.values()))
    composite_threshold = threshold_otsu(composite)
    return pd.DataFrame(
        [
            {
                "channel": channel,
                "mean": np.mean(image),
                "otsu_mean": np.mean(image[threshold_otsu(image) <= image]),
                "median": np.median(image),
                "otsu_median": np.median(image[threshold_otsu(image) <= image]),
                "composite_otsu_mean": np.mean(image[composite_threshold <= composite]),
                "composite_otsu_median": np.median(
                    image[composite_threshold <= composite]
                ),
                "composite_otsu_p90": np.percentile(
                    image[composite_threshold <= composite], 90
                ),
                "composite_otsu_p95": np.percentile(
                    image[composite_threshold <= composite], 95
                ),
                "composite_otsu_p98": np.percentile(
                    image[composite_threshold <= composite], 98
                ),
                "composite_otsu_p99": np.percentile(
                    image[composite_threshold <= composite], 99
                ),
                "p90": np.percentile(image, 90),
                "p95": np.percentile(image, 95),
                "p98": np.percentile(image, 98),
                "p99": np.percentile(image, 99),
                "p99.5": np.percentile(image, 99.5),
            }
            for channel, image in images.items()
        ],
    )


class Pipeline:
    # TODO: move boilerplate here

    def handle_message(self, msg):
        self._handle_message(msg)
        self._queue.poll()


class DefaultPipeline(Pipeline):
    DEFAULT_CONFIG = {
        "roi_detection_composite_func": power_law_composite,
        "segmentation_composite_func": partial(np.nanmean, axis=0),
        "roi_detection_func": find_trenches,
        "correct_drift": True,
        "segmentation_func": watershed_segment,
        "segment": True,
        "measure": True,
        "measure_fish": True,
    }
    REQUIRED_CONFIG = ["segmentation_channels"]
    FOV_T_SCHEMA = [("fov_num", pa.uint16()), ("t", pa.uint16())]
    FOV_T_ROI_SCHEMA = [
        ("fov_num", pa.uint16()),
        ("t", pa.uint16()),
        ("roi", pa.uint16()),
    ]
    FOV_CHANNEL_T_ROI_SCHEMA = [
        ("fov_num", pa.uint16()),
        ("channel", pa.string()),
        ("t", pa.uint16()),
        ("roi", pa.uint16()),
    ]

    def __init__(self, output_dir, config=None, delayed=False):
        self.output_dir = Path(output_dir)
        if not config:
            config = {}
        self.config = config
        self._apply_default_config()
        self._queue = DelayedQueue(wrapper=get_delayed(delayed))
        # fov
        self.first_t = {}
        # fov
        self.last_t = {}
        # fov
        self.fish_first_t = {}
        # fov, t
        self.initial_rois = DelayedTableStore(
            self._queue,
            output_dir / "initial_rois",
            schema=self.FOV_T_SCHEMA,
            write_options=dict(partition_cols=["fov_num", "t"]),
        )
        # fov, t
        self.rois = DelayedTableStore(
            self._queue,
            output_dir / "rois",
            schema=self.FOV_T_SCHEMA,
            write_options=dict(partition_cols=["fov_num", "t"]),
        )
        # fov, t
        self.fish_rois = DelayedTableStore(
            self._queue,
            output_dir / "fish_rois",
            schema=self.FOV_T_SCHEMA,
            write_options=dict(partition_cols=["fov_num", "t"]),
        )
        # fov, channel, t
        self.measurements = DelayedTableStore(
            self._queue,
            output_dir / "measurements",
            schema=self.FOV_CHANNEL_T_ROI_SCHEMA,
            write_options=dict(partition_cols=["fov_num", "t", "channel"]),
        )
        # fov, t
        self.mask_measurements = DelayedTableStore(
            self._queue,
            output_dir / "mask_measurements",
            schema=self.FOV_T_ROI_SCHEMA,
            write_options=dict(partition_cols=["fov_num", "t"]),
        )
        # fov, channel, t
        self.raw_frames = DelayedBatchedZarrStore(
            self._queue,
            output_dir / "raw_frames/fov={}",
            write_options=dict(chunks=(1, 1)),
        )
        # fov, channel, t
        self.processed_frames = DelayedBatchedZarrStore(
            self._queue,
            output_dir / "processed_frames/fov={}",
            write_options=dict(chunks=(1, 1)),
        )
        # fov, channel, t; roi
        self.crops = DelayedBatchedZarrStore(
            self._queue,
            output_dir / "crops/fov={}/channel={}",
            write_options=dict(chunks=(20, 40), incomplete_chunks_axes=[False, True]),
        )
        # fov, t; roi
        self.segmentation_crops = DelayedBatchedZarrStore(
            self._queue,
            output_dir / "segmentation_crops/fov={}",
            write_options=dict(chunks=(20, 40), incomplete_chunks_axes=[False, True]),
        )
        # fov, t; roi
        self.segmentation_masks = DelayedBatchedZarrStore(
            self._queue,
            output_dir / "segmentation_masks/fov={}",
            write_options=dict(chunks=(20, 40), incomplete_chunks_axes=[False, True]),
        )
        # fov
        self.initial_drift_features = DelayedStore(self._queue)
        # fov
        self.image_limits = DelayedStore(self._queue)
        # fov
        self.fish_image_limits = DelayedStore(self._queue)
        # fov, t
        self.shifts = DelayedStore(self._queue)
        # fov, t
        self.fish_shifts = DelayedStore(self._queue)
        # fov, channel, t
        self.fish_raw_frames = DelayedBatchedZarrStore(
            self._queue,
            output_dir / "fish_raw_frames/fov={}",
            write_options=dict(chunks=(1, 1)),
        )
        # fov, channel, t
        self.fish_processed_frames = DelayedBatchedZarrStore(
            self._queue,
            output_dir / "fish_processed_frames/fov={}",
            write_options=dict(chunks=(1, 1)),
        )
        # fov, channel, t; roi
        self.fish_crops = DelayedBatchedZarrStore(
            self._queue,
            output_dir / "fish_crops/fov={}/channel={}",
            write_options=dict(chunks=(10, 40), incomplete_chunks_axes=[False, True]),
        )
        # fov, t
        self.fish_measurements = DelayedTableStore(
            self._queue,
            output_dir / "fish_measurements",
            schema=self.FOV_T_ROI_SCHEMA,
            write_options=dict(partition_cols=["fov_num", "t"]),
        )
        self._stores = [x for x in vars(self).values() if isinstance(x, DelayedStore)]

    def _apply_default_config(self):
        config = {**self.DEFAULT_CONFIG, **self.config}
        config["trench_detection_channels"] = config.get(
            "trench_detection_channels"
        ) or config.get("segmentation_channels")
        if config["measure"] and not config["segment"]:
            raise ValueError("cannot enable measurement with segmentation disabled")
        self.config = config

    def delayed(self, func, *args, **kwargs):
        return self._queue.delayed(func, *args, **kwargs)

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
                self.handle_done()
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
        # print("IMAGE", metadata)
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
        roi_detection_frame = self.delayed(
            self.config["roi_detection_composite_func"], segmentation_frames
        )
        # segmentation_frame = self.delayed(
        #     self.config["segmentation_composite_func"], segmentation_frames
        # )
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
                    roi_detection_frame,
                    _keep_args=True,
                    **(self.config.get("roi_detection_kwargs") or {}),
                ),
            )
            shift = self.shifts.setdefault((fov_num, t), np.array([0, 0]))
            self.initial_drift_features.setdefault(
                fov_num,
                self.delayed(get_drift_features, roi_detection_frame, initial_rois),
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
                    roi_detection_frame,
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
        # segmentation_crops
        segmentation_crops = self.segmentation_crops.setdefault(
            (fov_num, t),
            self.delayed(
                combine_crops,
                self.config["segmentation_composite_func"],
                [
                    self.crops[(fov_num, channel, t)]
                    for channel in self.config["segmentation_channels"]
                ],
            ),
        )
        # if crop available -> segment crops
        if self.config["segment"]:
            seg_masks = self.segmentation_masks.setdefault(
                (fov_num, t),
                self.delayed(
                    segment_crops,
                    segmentation_crops,
                ),
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
        # cleanup
        self.write_stores()
        for store in [
            self.raw_frames,
            self.processed_frames,
            self.crops,
            self.segmentation_crops,
            self.segmentation_masks,
            self.fish_raw_frames,
            self.fish_processed_frames,
            self.fish_crops,
            self.measurements,
            self.mask_measurements,
            self.fish_measurements,
            self.rois,
        ]:
            for key in list(store.keys()):
                if key[-1] != t:
                    del store[key]

    def handle_fish_barcode(self, msg):
        image = msg["image"]
        metadata = msg["metadata"]
        # print("FISH IMAGE", metadata)
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
        fish_image_limits = self.fish_image_limits.setdefault(
            fov_num,
            self.delayed(lambda img: get_image_limits(img.shape), processed_frame),
        )
        fish_first_t = self.fish_first_t.setdefault(fov_num, t)
        last_t = self.last_t.get(fov_num, None)
        if last_t is None:
            raise ValueError("received FISH image before ROI detection")
        last_science_rois = self.rois[(fov_num, last_t)]
        if t == fish_first_t:
            rois = self.fish_rois.setdefault((fov_num, t), last_science_rois)
        else:
            drift_tracking_channel = self.config["fish_drift_tracking_channel"]
            shift = self.fish_shifts.setdefault(
                (fov_num, t),
                self.delayed(
                    # phase_cross_correlation returns offset in y,x (array) order
                    # not x,y order as filter_and_shift_rois expects
                    # also we swap ref and moving image (equivalently, we could've inverted sign)
                    # to get correct offset sign
                    compose(
                        lambda x: x[0].astype(np.int16)[::-1], phase_cross_correlation
                    ),
                    processed_frames[(fov_num, drift_tracking_channel, t)],
                    processed_frames[(fov_num, drift_tracking_channel, fish_first_t)],
                ),
            )
            rois = self.fish_rois.setdefault(
                (fov_num, t),
                self.delayed(
                    filter_and_shift_rois,
                    last_science_rois,
                    shift,
                    fish_image_limits,
                ),
            )
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
        # cleanup
        # flush non-FISH DelayedBatchedZarrStores so that non-FISH write tasks launch before FISH tasks
        self.crops.flush()
        self.segmentation_crops.flush()
        self.segmentation_masks.flush()
        self.write_stores()
        if self.config.get("preprocess_func"):
            for key in list(self.fish_raw_frames):
                del self.fish_raw_frames[key]
            for key in list(self.fish_processed_frames):
                if key[0] == fov_num and key[1:3] != (
                    self.config["fish_drift_tracking_channel"],
                    fish_first_t,
                ):
                    del self.fish_processed_frames[key]
        else:
            for key in list(self.fish_raw_frames):
                if key[0] == fov_num and key[1:3] != (
                    self.config["fish_drift_tracking_channel"],
                    fish_first_t,
                ):
                    del self.fish_raw_frames[key]
        for store in [self.fish_crops, self.fish_measurements]:
            for key in list(store.keys()):
                if key[0] == fov_num and key[-1] != t:
                    del store[key]
        for store in [
            self.raw_frames,
            self.processed_frames,
            self.crops,
            self.segmentation_crops,
            self.segmentation_masks,
            self.measurements,
            self.mask_measurements,
        ]:
            store.flush()
            for key in list(store.keys()):
                del store[key]
        for store in [self.rois]:
            for key in list(store.keys()):
                if key[0] == fov_num and key[1] != last_t:
                    del store[key]
        for store in [self.fish_rois]:
            for key in list(store.keys()):
                if key[0] == fov_num and key[1] != t:
                    del store[key]

    def write_stores(self):
        ####
        # for store in self._stores:
        #     store.write()
        #### DelayedTableStore
        self.initial_rois.write()
        self.rois.write()
        self.measurements.write()
        self.mask_measurements.write()
        self.fish_measurements.write()
        #### DelayedBatchedZarrStore
        # don't write, but allow entries to be deleted
        for store in [
            self.raw_frames,
            self.processed_frames,
            self.fish_raw_frames,
            self.fish_processed_frames,
        ]:
            store._write_queue.clear()
        # self.raw_frames.write()
        # self.processed_frames.write()
        self.crops.write()
        self.segmentation_crops.write()
        self.segmentation_masks.write()
        # self.fish_raw_frames.write()
        # self.fish_processed_frames.write()
        self.fish_crops.write()

    def handle_done(self):
        print("CLEANING UP")
        # for store in self._stores:
        #     store.write()
        # self.write_stores() # we should already have finished writing after every handle_message
        for store in self._stores:
            store.flush()
            for key in list(store.keys()):
                del store[key]
        print("DONE")
