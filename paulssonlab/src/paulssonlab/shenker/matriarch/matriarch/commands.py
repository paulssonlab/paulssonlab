import pandas as pd
from distributed import Client
from dask_jobqueue import SLURMCluster
import click

# from .processing import compress_trenches
from .util import multi_join, apply_map_futures

IDX = pd.IndexSlice


@click.group()
def dask_command():
    pass


@dask_command.command()
def compress_nd2(in_path):
    all_frames, metadata, parsed_metadata = workflow.get_nd2_frame_list(nd2_filenames)
    image_limits = workflow.get_filename_image_limits(metadata)
    frames_to_process = all_frames.loc[IDX[:, :, ["MCHERRY"], 0], :]
    find_trenches_diag = diagnostics.wrap_diagnostics(
        trench_detection.hough.find_trenches, ignore_exceptions=True, pandas=True
    )
    trench_info_futures = {
        idx: client.submit(
            find_trenches_diag, client.submit(workflow.get_nd2_frame, **idx._asdict())
        )
        for idx, row in util.iter_index(frames_to_process)
    }
    cluster = SLURMCluster(
        queue="short",
        walltime="03:00:00",
        # job_extra=['--qos=HHH'],
        memory="1GB",
        local_directory="/tmp",
        cores=1,
        processes=1,
        env_extra=[
            'export PYTHONPATH="/home/jqs1/projects/matriarch"',
            'export LD_PRELOAD="/home/jqs1/lib/libjemalloc.so.2"',
        ],
    )
    client = Client(cluster)
    trench_info = util.apply_map_futures(
        client.gather, trench_info_futures, predicate=lambda x: x.status == "finished"
    )
    trench_points, trench_diag, trench_err = workflow.unzip_trench_info(trench_info)
    bad_pitch = (
        trench_diag["find_trench_lines.hough_2.peak_func.pitch"] - 24
    ).abs() > 1
    trench_points_good = trench_points[~multi_join(trench_points.index, bad_pitch)]
    trench_bbox_futures = []
    for _, trenches in trench_points_good.groupby(["filename", "position", "t"]):
        trench_bbox_futures.append(
            client.submit(workflow.get_trench_bboxes, trenches, image_limits)
        )
        trench_bbox_results = util.apply_map_futures(
            client.gather,
            trench_bbox_futures,
            predicate=lambda x: x.status == "finished",
        )
        trench_bboxes = pd.concat(
            [trench_points_good, pd.concat(trench_bbox_results, axis=0)], axis=1
        )
        trench_bboxes_t0 = util.get_one(trench_bboxes.groupby("t"))[1]
        selected_trenches_segmentation = trench_bboxes_t0[
            trench_bboxes_t0[("info", "hough_value")] > 90
        ].loc[IDX[:, :, ["MCHERRY"], 0, :, :], :]
        selected_trenches_segmentation.index = (
            selected_trenches_segmentation.index.droplevel("channel")
        )
        frames_to_analyze = all_frames.loc[IDX[:, :, ["MCHERRY", "YFP"], :], :]
        zarrify = partial(
            zarr.array,
            compressor=processing.DEFAULT_COMPRESSOR,
            order=processing.DEFAULT_ORDER,
            chunks=False,
        )
        frame_transformation = compose(
            zarrify,
            partial(image.quantize, bits=8),
            partial(image.downsample, factor=4),
        )
        compress_frames(
            selected_trenches_segmentation2,
            frames_to_analyze2,
            include_frame=True,
            frame_transformation=frame_transformation,
        )


commands = [compress_nd2]
