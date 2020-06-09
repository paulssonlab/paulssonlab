#!/usr/bin/env python3

import click

from convert import main_conversion_function
from trench_detect_new import main_detection_function
from segmentation import main_segmentation_function

# Modified from https://rosettacode.org/wiki/Range_expansion#Python
# Input a range, such as: 1,4-7
# Return a list with all elements within the range.
def range_expand(range_string):
    result = []

    for r in range_string.split(","):
        # Start at 1, because the first number might be negative,
        # meaning that the first char is a minus (dash) sign.
        if "-" in r[1:]:
            first, second = r[1:].split("-", 1)
            result += range(int(r[0] + first), int(second) + 1)

        else:
            result.append(int(r))

    return result


# Make an empty cli group which doesn't do anything, but which hosts the sub-commands.
@click.group()
def cli():
    pass


@cli.command(help="Convert a directory of ND2 files to an HDF5 file.")
@click.option(
    "-o",
    "--out-dir",
    "out_dir",
    required=True,
    default="HDF5",
    type=str,
    help="Output directory.",
    show_default=True,
)
@click.option(
    "-i",
    "--in-dir",
    "in_dir",
    required=True,
    default="ND2",
    type=str,
    help="Input directory of ND2 files.",
    show_default=True,
)
@click.option(
    "-n",
    "--num-cpu",
    "num_cpu",
    required=False,
    default=None,
    type=click.IntRange(1, None, clamp=True),
    help="Number of CPUs to use. [default: all CPUs]",
    show_default=False,
)
@click.option(
    "-f",
    "--frames",
    "frames",
    required=False,
    default=None,
    type=str,
    help="List of frames to copy. [default: all frames]",
    show_default=False,
)
@click.option(
    "-F",
    "--fovs",
    "fovs",
    required=False,
    default=None,
    type=str,
    help="List of FOVs to copy. [default: all FOVs]",
    show_default=False,
)
def convert(out_dir, in_dir, num_cpu, frames, fovs):

    # a. Range of time frames to be analyzed, if the user doesn't want to copy all the frames
    if frames:
        frames = range_expand(frames)

    # b. Range of FOVs to be analyzed
    if fovs:
        fovs = range_expand(fovs)

    main_conversion_function(out_dir, in_dir, num_cpu, frames, fovs)


@cli.command(
    help="Detect trenches and write their rectangular regions to an HDF5  file."
)
@click.option(
    "-o",
    "--out-dir",
    "out_dir",
    required=True,
    default="REGIONS",
    type=str,
    help="Output directory.",
    show_default=True,
)
@click.option(
    "-i",
    "--in-file",
    "in_file",
    required=True,
    default="HDF5/data.h5",
    type=str,
    help="Input HDF5 file.",
    show_default=True,
)
@click.option(
    "-n",
    "--num-cpu",
    "num_cpu",
    required=False,
    default=None,
    type=click.IntRange(1, None, clamp=True),
    help="Number of CPUs to use. [default: all CPUs]",
    show_default=False,
)
@click.option(
    "-P",
    "--params-file",
    "params_file",
    required=True,
    default="seg_params.yaml",
    type=str,
    help="Regions detection parameters file (YAML).",
    show_default=True,
)
@click.option(
    "s",
    "--share-regions",
    required=True,
    default=False,
    type=bool,
    help="Share region detection across frames (detect only within the first frame)",
)
def trench_detect(out_dir):
    main_detection_function(out_dir, in_file, num_cpu, params_file, share_regions)


@cli.command(help="Detect cells and write their properties and masks to an HDF5 file.")
@click.option(
    "-o",
    "--out-dir",
    "out_dir",
    required=True,
    default="SEGMENTATION",
    type=str,
    help="Output for new HDF5 directory for masks & measurements.",
    show_default=True,
)
@click.option(
    "-i",
    "--in-filer",
    "in_file",
    required=True,
    default="HDF5/data.h5",
    type=str,
    help="Input HDF5 file with images.",
    show_default=True,
)
@click.option(
    "-n",
    "--num-cpu",
    "num_cpu",
    required=False,
    default=None,
    type=click.IntRange(1, None, clamp=True),
    help="Number of CPUs to use. [default: all CPUs]",
    show_default=False,
)
@click.option(
    "-P",
    "--params-file",
    "params_file",
    required=True,
    default="seg_params.yaml",
    type=str,
    help="Segmentation parameters file (YAML).",
    show_default=True,
)
@click.option(
    "-R",
    "--regions-file",
    "regions_file",
    required=False,
    default=None,
    type=str,
    help="HDF5 file containing image regions. [default: analyze entire image, no regions]",
    show_default=True,
)  # FIXME what should the default be???
def segment(out_dir, in_file, num_cpu, params_file, regions_file):
    main_segmentation_function(out_dir, in_file, num_cpu, params_file, regions_file)


###

if __name__ == "__main__":
    cli()
