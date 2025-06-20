#!/usr/bin/env python3

import click
from convert import main_conversion_function
from corrections import main_corrections_function
from lineages import main_lineages_function
from metadata import main_print_metadata_function
from napari_browse_hdf5 import main_hdf5_browser_function
from napari_browse_kymographs import main_kymograph_browser_function
from napari_browse_nd2 import main_nd2_browser_function
from renumber_trenches import main_renumbering_function
from segmentation import main_segmentation_function
from trench_detect import main_detection_function
from trench_measurements import main_trench_measurements_function
from write_kymographs import main_kymographs_function


def range_expand(range_string):
    """
    Input a range, such as: 1,4-7

    Return a list with all elements within the range.

    Modified from https://rosettacode.org/wiki/Range_expansion#Python
    """
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


@click.group()
def cli():
    """Invoke a Command to perform the desired operation:"""
    pass


@cli.command(no_args_is_help=True)
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
    "-S",
    "--settings-file",
    "settings_file",
    required=True,
    default="trenchcoat_settings.yaml",
    type=str,
    help="TrenchCoat settings file (YAML).",
    show_default=True,
)
def browse_nd2(in_dir, settings_file):
    """Use Napari to browse a directory of ND2 files."""
    main_nd2_browser_function(in_dir=in_dir, settings_file=settings_file)


@cli.command(no_args_is_help=True)
@click.option(
    "-i",
    "--images-file",
    "images_file",
    required=True,
    default="HDF5/",
    type=str,
    help="Input HDF5 file with images.",
    show_default=True,
)
@click.option(
    "-m",
    "--masks-file",
    "masks_file",
    required=False,
    type=str,
    help="Input HDF5 file with masks.",
)
@click.option(
    "-r",
    "--regions-file",
    "regions_file",
    required=False,
    type=str,
    help="Input HDF5 file with regions.",
)
@click.option(
    "-S",
    "--settings-file",
    "settings_file",
    required=True,
    default="trenchcoat_settings.yaml",
    type=str,
    help="TrenchCoat settings file (YAML).",
    show_default=True,
)
@click.option(
    "-V",
    "--viewer-params",
    "viewer_params_file",
    required=False,
    type=str,
    help="YAML file with image viewer params & information for displaying column data.",
)
@click.option(
    "-d",
    "--data-table-file",
    "data_table_file",
    required=False,
    type=str,
    help="Path to HDF5 file containing cell measurements for computed image.",
)
def browse_hdf5(
    images_file,
    masks_file,
    regions_file,
    settings_file,
    viewer_params_file,
    data_table_file,
):
    """
    Use Napari to browse a dataset & to visualize trenches and cell masks.

    Add option to "compute" an image based on properties within segmented regions.
    """
    main_hdf5_browser_function(
        images_file=images_file,
        masks_file=masks_file,
        regions_file=regions_file,
        settings_file=settings_file,
        viewer_params_file=viewer_params_file,
        data_table_file=data_table_file,
    )


@cli.command(no_args_is_help=True)
@click.option(
    "-i",
    "--images-file",
    "images_file",
    required=True,
    default="HDF5/",
    type=str,
    help="Input HDF5 file with images.",
    show_default=True,
)
@click.option(
    "-m",
    "--masks-file",
    "masks_file",
    required=False,
    type=str,
    help="Input HDF5 file with masks.",
)
@click.option(
    "-r",
    "--regions-file",
    "regions_file",
    required=True,
    type=str,
    default="regions_table.h5",
    help="Input HDF5 file with regions (trenches must be re-labeled for left/right drift).",
)
@click.option(
    "-L",
    "--lineages-file",
    "lineages_file",
    required=False,
    type=str,
    help="Input HDF5 file with cell measurements, re-labeled trenches and cell lineages.",
)
@click.option(
    "-S",
    "--settings-file",
    "settings_file",
    required=True,
    default="trenchcoat_settings.yaml",
    type=str,
    help="TrenchCoat settings file (YAML).",
    show_default=True,
)
@click.option(
    "-V",
    "--viewer-params",
    "viewer_params_file",
    required=False,
    type=str,
    help="YAML file with image viewer params & information for displaying column data.",
)
@click.option(
    "-d",
    "--data-table-file",
    "data_table_file",
    required=False,
    type=str,
    help="Path to HDF5 file containing cell measurements for computed image.",
)
@click.option(
    "-P",
    "--precomputed",
    "precomputed",
    required=True,
    default=False,
    type=bool,
    help="Specify whether kymographs of intensity images & masks were pre-computed.",
)
def browse_kymographs(
    images_file,
    masks_file,
    regions_file,
    settings_file,
    lineages_file,
    viewer_params_file,
    data_table_file,
    precomputed,
):
    """
    Use Napari to browse kymographs.
    The regions file must be post-processed for correcting stage drift,
    and for merging cell information with trench information.
    Corrections are unimplemented.
    """
    main_kymograph_browser_function(
        images_file=images_file,
        masks_file=masks_file,
        regions_file=regions_file,
        settings_file=settings_file,
        lineages_file=lineages_file,
        viewer_params_file=viewer_params_file,
        data_table_file=data_table_file,
        precomputed=precomputed,  # Are the intensity images, masks pre-computed?
    )


@cli.command(no_args_is_help=True)
@click.option(
    "-o",
    "--out-dir",
    "out_dir",
    required=True,
    default="HDF5/",
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
    help="Number of CPUs to use.  [default: all CPUs]",
    show_default=False,
)
@click.option(
    "-f",
    "--frames",
    "frames",
    required=False,
    default=None,
    type=str,
    help="List of frames to copy.  [default: all frames]",
    show_default=False,
)
@click.option(
    "-F",
    "--fovs",
    "fovs",
    required=False,
    default=None,
    type=str,
    help="List of FOVs to copy.  [default: all FOVs]",
    show_default=False,
)
def convert(out_dir, in_dir, num_cpu, frames, fovs):
    """Convert a directory of ND2 files to an HDF5 file."""
    # a. Range of time frames to be analyzed, if the user doesn't want to copy all the frames
    if frames:
        frames = range_expand(frames)

    # b. Range of FOVs to be analyzed
    if fovs:
        fovs = range_expand(fovs)

    main_conversion_function(
        hdf_dir=out_dir,
        nd2_dir=in_dir,
        num_cpu=num_cpu,
        frames=frames,
        fields_of_view=fovs,
    )


@cli.command(no_args_is_help=True)
@click.option(
    "-o",
    "--out-file",
    "out_file",
    required=True,
    default="regions_table.h5",
    type=str,
    help="Output directory.",
    show_default=True,
)
@click.option(
    "-i",
    "--in-file",
    "in_file",
    required=True,
    default="HDF5/",
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
    default="trench_params.yaml",
    type=str,
    help="Regions detection parameters file (YAML).",
    show_default=True,
)
def trench_detect(out_file, in_file, num_cpu, params_file):
    """Detect trenches and write their rectangular regions to an HDF5 file."""
    main_detection_function(
        out_file=out_file, in_file=in_file, num_cpu=num_cpu, params_file=params_file
    )


@cli.command(no_args_is_help=True)
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
    "--in-file",
    "in_file",
    required=True,
    default="HDF5/",
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
    """Detect cells and write their properties and masks to an HDF5 file."""
    main_segmentation_function(
        out_dir=out_dir,
        in_file=in_file,
        num_cpu=num_cpu,
        params_file=params_file,
        regions_file=regions_file,
    )


@cli.command(no_args_is_help=True)
@click.option(
    "-o",
    "--out-dir",
    "out_dir",
    required=True,
    default="TRENCH_MEASUREMENTS",
    type=str,
    help="Output for new HDF5 directory for trench measurements.",
    show_default=True,
)
@click.option(
    "-i",
    "--in-file",
    "in_file",
    required=True,
    default="HDF5/",
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
    "-R",
    "--regions-file",
    "regions_file",
    required=False,
    default=None,
    type=str,
    help="HDF5 file containing image regions. [default: analyze entire image, no regions]",
    show_default=True,
)  # FIXME what should the default be???
def trench_measurements():
    """
    Analyze whole trenches, without cell segmentation.
    FIXME This has not been tested!
    """
    main_trench_measurements_function(
        out_dir=out_dir, in_file=in_file, num_cpu=num_cpu, regions_file=regions_file
    )


@cli.command(no_args_is_help=True)
@click.option(
    "-o",
    "--out-file",
    "out_file",
    required=True,
    default="corrections.h5",
    type=str,
    help="Output HDF5 file with corrections matrices.",
    show_default=True,
)
@click.option(
    "-i",
    "--in-file",
    "in_file",
    required=True,
    default="BLANKS/data.h5",
    type=str,
    help="Input HDF5 file with images.",
    show_default=True,
)
@click.option(
    "-D",
    "--dark-channel",
    "dark_channel",
    required=False,
    default=None,
    type=str,
    help="Name of dark channel (cannot be used with -B).",
    show_default=False,
)
@click.option(
    "-B",
    "--background-values-file",
    "bg_file",
    required=False,
    default=None,
    type=str,
    help="YAML file with background values (cannot be used with -D).",
    show_default=False,
)
def corrections(in_file, out_file, dark_channel, bg_file):
    """Generate camera bias and flat field corrections matrices from images."""
    main_corrections_function(
        in_file=in_file, out_file=out_file, dark_channel=dark_channel, bg_file=bg_file
    )


@cli.command(no_args_is_help=True)
@click.option(
    "-i",
    "--in-file",
    "in_file",
    required=True,
    default="HDF5/data.h5",
    type=str,
    help="Input HDF5 or ND2 file with metadata.",
    show_default=True,
)
@click.option(
    "-f",
    "--sub-file-name",
    "sub_file_name",
    required=False,
    default=None,
    type=str,
    help="Name of sub-file to print metadata for. This file was a separate ND2 file before conversion to HDF5.",
    show_default=False,
)
def print_metadata(in_file, sub_file_name):
    """Print the HDF5 metadata to the terminal. If no sub-file is specified, then print a list of all sub-files."""
    main_print_metadata_function(in_file=in_file, sub_file_name=sub_file_name)


@cli.command(no_args_is_help=True)
@click.option(
    "-i",
    "--in-file",
    "in_file",
    required=True,
    default="HDF5/",
    type=str,
    help="Input HDF5 file with metadata.",
    show_default=True,
)
@click.option(
    "-R",
    "--regions-file",
    "regions_file",
    required=False,
    default="regions_table.h5",
    type=str,
    help="HDF5 file containing image regions",
    show_default=True,
)
@click.option(
    "-O",
    "--regions-out-file",
    "regions_outfile",
    required=True,
    default=None,
    type=str,
    help="Write renumbered trenches & regions.",
    show_default=True,
)
@click.option(
    "-M",
    "--measurements-file",
    "cell_measurements_file",
    required=False,
    default="SEG/",
    type=str,
    help="Directory containing subdir with cell measurements with original trench numbering.",
    show_default=True,
)
@click.option(
    "-m",
    "--measurements-out-file",
    "cell_measurements_outfile",
    required=False,
    default=None,
    type=str,
    help="Write renumbered trenches paired with cell measurements.",
    show_default=True,
)
@click.option(
    "-H",
    "--method",
    "method",
    required=True,
    default="stage",
    type=str,
    help="Method to use for calculating drift. Options are stage or image.",
    show_default=True,
)
@click.option(
    "-C",
    "--channel",
    "channel",
    required=False,
    default=None,
    type=str,
    help="Channel to use if using images to correct for drift.",
    show_default=True,
)
# NOTE The value of this threshold working depends on how much drift there is left unaccounted for
# by the stage offsets.
# - A lower value be more likely to cause re-numbering.
# - An intermediate value could have the weird effect of only re-numbering some
# of the frames.
# - A high value will not have any effect.
@click.option(
    "-T",
    "--threshold",
    "threshold",
    required=True,
    default=10,
    type=int,
    help="Threshold for calling a renumbering event (pixels).",
    show_default=True,
)
def correct_trench_numbering(
    in_file,
    regions_file,
    regions_outfile,
    cell_measurements_file,
    cell_measurements_outfile,
    threshold,
    method,
    channel,
):
    """Take into account when a trench lies near the edge of an image by re-numbering all the trenches as needed."""
    main_renumbering_function(
        in_file=in_file,
        regions_file=regions_file,
        regions_outfile=regions_outfile,
        cell_measurements_file=cell_measurements_file,
        cell_measurements_outfile=cell_measurements_outfile,
        threshold=threshold,
        method=method,
        channel=channel,
    )


@cli.command(no_args_is_help=True)
@click.option(
    "-i",
    "--in-file",
    "in_file",
    required=True,
    default="newsegs.h5",  # TODO come up with a better name
    type=str,
    help="Input HDF5 file with cell measurements data and trench labels.",
    show_default=True,
)
@click.option(
    "-o",
    "--out-file",
    "out_file",
    required=True,
    default="lineages.h5",
    type=str,
    help="Output HDF5 file with lineages table.",
    show_default=True,
)
@click.option(
    "-L",
    "--length-buffer",
    "length_buffer",
    required=True,
    default=5,
    type=int,
    help="Fudge factor: a change in length less than this much is not considered to reflect a cell division.",
    show_default=True,
)
@click.option(
    "-T",
    "--trench-length",
    "trench_length",
    required=True,
    type=str,
    help="Length of a trench. Used to check whether a cell is at the trench exit.",
)
def lineage_tracking(in_file, out_file, length_buffer, trench_length):
    """Track lineages using HDF5 tables. Write a new HDF5 table."""
    main_lineages_function(
        infile=in_file,
        outfile=out_file,
        length_buffer=length_buffer,
        trench_length=trench_length,
    )


@cli.command(no_args_is_help=True)
@click.option(
    "-i",
    "--images-file",
    "images_file",
    required=True,
    default="HDF5/data.h5",
    type=str,
    help="Input HDF5 file with intensity images. Also required for converting masks, because it contains essential metadata.",
    show_default=True,
)
@click.option(
    "-m",
    "--masks-file",
    "masks_file",
    required=False,
    type=str,
    help="Input HDF5 file with segmentation masks.",
    show_default=True,
)
@click.option(
    "-R",
    "--regions-file",
    "regions_file",
    required=True,
    default="corrected_regions.h5",
    type=str,
    help='HDF5 file with regions coordinates table. Assumes that stage drift has been accounted for and that trenches have been renumbered ("corrected_trench_label").',
    show_default=True,
)
@click.option(
    "-o",
    "--out-dir",
    "out_dir",
    required=True,
    type=str,
    help="Directory to write HDF5 file with kymographs.",
)
@click.option(
    "-n",
    "--num-cpu",
    "num_cpu",
    required=False,
    default=None,
    type=click.IntRange(1, None, clamp=True),
    help="Number of CPUs to use.  [default: all CPUs]",
    show_default=False,
)
def write_kymographs(images_file, masks_file, regions_file, out_dir, num_cpu):
    """Pre-render intensity image kymographs & write to disk."""
    main_kymographs_function(
        images_file=images_file,
        masks_file=masks_file,
        regions_file=regions_file,
        out_dir=out_dir,
        num_cpu=num_cpu,
    )


###

if __name__ == "__main__":
    cli()
